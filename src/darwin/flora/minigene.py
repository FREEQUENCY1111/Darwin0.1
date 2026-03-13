"""
MiniGeneHunter — Small gene recovery organism.

Feeds on: genes.called (to know what's already predicted)
Produces: minigenes.found (small ORFs with RBS motifs)

Like mycorrhizal fungi that find nutrients in places roots
can't reach — this organism finds tiny genes that Prodigal
misses (ORFs < 150bp with strong Shine-Dalgarno signals).
"""

from __future__ import annotations

import logging

from darwin.flora.base import Organism
from darwin.rocks.models import Feature, FeatureType, Genome, Strand
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.flora.minigene")

# Shine-Dalgarno consensus: AGGAGG (or subsets)
# Score each position in the motif
SD_MOTIFS = [
    ("AGGAGG", 6.0),
    ("AGGAG", 5.0),
    ("GGAGG", 5.0),
    ("AGGA", 4.0),
    ("GAGG", 4.0),
    ("GGA", 3.0),
    ("AGG", 3.0),
    ("GAG", 3.0),
]

# Standard start and stop codons
START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}

# Genetic code table 11 (bacterial)
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(comp)[::-1]


def _translate(nuc: str) -> str:
    """Translate nucleotide sequence to protein (table 11)."""
    protein = []
    for i in range(0, len(nuc) - 2, 3):
        codon = nuc[i : i + 3].upper()
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def _score_rbs(upstream: str) -> float:
    """
    Score the Shine-Dalgarno motif in the upstream region.

    Scans 20bp upstream of start codon (positions -7 to -12
    being the sweet spot for SD interaction with 16S rRNA).
    """
    upstream = upstream.upper()
    best_score = 0.0
    for motif, score in SD_MOTIFS:
        if motif in upstream:
            best_score = max(best_score, score)
    return best_score


class MiniGeneHunter(Organism):
    """Small gene recovery — the mycorrhizal network."""

    name = "minigene_hunter"
    feeds_on_nutrients = [NutrientType.GENES_CALLED]
    produces_nutrients = [NutrientType.MINIGENES_FOUND]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return True  # Pure Python, no external tools needed

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """
        Scan intergenic regions for small ORFs with RBS motifs.

        Strategy: find ORFs 45-150bp in intergenic gaps, score
        their Shine-Dalgarno motif, keep those with a strong SD
        signal (score >= 5.0) or a moderate SD + ATG start codon.
        This balances sensitivity for real leader peptides against
        precision by filtering chance 3-mer motif hits.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        locus_prefix = config.get("locus_tag_prefix", "DARWIN")

        self.logger.info("🔍 Hunting for small genes in intergenic regions...")

        # Get existing features to find intergenic gaps
        existing_features = genome.all_features
        existing_tag_nums = set()
        for f in existing_features:
            if f.locus_tag and f.locus_tag.startswith(locus_prefix + "_"):
                try:
                    num = int(f.locus_tag.split("_")[-1])
                    existing_tag_nums.add(num)
                except ValueError:
                    pass

        # Start numbering after the highest existing tag
        next_num = max(existing_tag_nums, default=0) + 1
        minigene_count = 0
        min_orf_bp = 45         # raised from 30 — <15 AA ORFs are mostly noise
        max_orf_bp = 150
        min_rbs_score = 5.0     # raised from 4.0 — require 5-mer SD motif minimum

        for contig in genome.contigs:
            seq = contig.sequence.upper()
            if len(seq) < min_orf_bp:
                continue

            # Build occupied intervals from existing features
            occupied = []
            for f in contig.features:
                occupied.append((f.start, f.end))
            occupied.sort()

            # Find intergenic gaps
            gaps = self._find_gaps(occupied, len(seq))

            # Scan each gap on both strands
            for gap_start, gap_end in gaps:
                gap_len = gap_end - gap_start
                if gap_len < min_orf_bp:
                    continue

                for strand in [Strand.FORWARD, Strand.REVERSE]:
                    if strand == Strand.FORWARD:
                        scan_seq = seq[gap_start - 1 : gap_end]
                    else:
                        scan_seq = _reverse_complement(seq[gap_start - 1 : gap_end])

                    orfs = self._find_orfs(scan_seq, min_orf_bp, max_orf_bp)

                    for orf_start_rel, orf_end_rel in orfs:
                        # Get upstream region for RBS scoring
                        upstream_start = max(0, orf_start_rel - 20)
                        upstream = scan_seq[upstream_start:orf_start_rel]
                        rbs_score = _score_rbs(upstream)

                        if rbs_score < min_rbs_score:
                            continue

                        # Non-ATG starts (GTG/TTG) in small ORFs are mostly
                        # noise — require the strongest SD signal (6-mer)
                        start_codon = scan_seq[orf_start_rel:orf_start_rel + 3].upper()
                        if start_codon != "ATG" and rbs_score < 6.0:
                            continue

                        # Convert to genomic coordinates
                        if strand == Strand.FORWARD:
                            feat_start = gap_start + orf_start_rel
                            feat_end = gap_start + orf_end_rel - 1
                        else:
                            feat_end = gap_end - orf_start_rel
                            feat_start = gap_end - orf_end_rel + 1

                        # Ensure within gap
                        if feat_start < gap_start or feat_end > gap_end:
                            continue

                        # Check no overlap with existing features
                        if self._overlaps_existing(feat_start, feat_end, occupied):
                            continue

                        # Translate
                        orf_nuc = scan_seq[orf_start_rel:orf_end_rel]
                        translation = _translate(orf_nuc)
                        if len(translation) < 10:  # at least 10 AA
                            continue

                        locus_tag = f"{locus_prefix}_{next_num:05d}"
                        next_num += 1

                        feature = Feature(
                            type=FeatureType.CDS,
                            start=feat_start,
                            end=feat_end,
                            strand=strand,
                            score=rbs_score,
                            contig_id=contig.id,
                            locus_tag=locus_tag,
                            product="hypothetical protein",
                            inference="ab initio prediction:Darwin:MiniGeneHunter",
                            translation=translation,
                            note=f"small_orf;rbs_score={rbs_score:.1f}",
                        )
                        contig.features.append(feature)
                        minigene_count += 1

        self.logger.info(f"🔬 Found {minigene_count} small genes with RBS motifs")

        return Nutrient(
            type=NutrientType.MINIGENES_FOUND,
            data={
                "genome": genome,
                "minigene_count": minigene_count,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    @staticmethod
    def _find_gaps(occupied: list[tuple[int, int]], seq_len: int) -> list[tuple[int, int]]:
        """Find intergenic gaps between occupied intervals."""
        if not occupied:
            return [(1, seq_len)]

        gaps = []
        # Gap before first feature
        if occupied[0][0] > 1:
            gaps.append((1, occupied[0][0] - 1))

        # Gaps between features
        for i in range(len(occupied) - 1):
            gap_start = occupied[i][1] + 1
            gap_end = occupied[i + 1][0] - 1
            if gap_end > gap_start:
                gaps.append((gap_start, gap_end))

        # Gap after last feature
        if occupied[-1][1] < seq_len:
            gaps.append((occupied[-1][1] + 1, seq_len))

        return gaps

    @staticmethod
    def _find_orfs(seq: str, min_len: int, max_len: int) -> list[tuple[int, int]]:
        """Find all ORFs (start to stop codon) within length range."""
        orfs = []
        seq = seq.upper()
        seq_len = len(seq)

        for frame in range(3):
            i = frame
            while i < seq_len - 2:
                codon = seq[i : i + 3]
                if codon in START_CODONS:
                    # Scan for stop codon
                    for j in range(i + 3, min(i + max_len + 3, seq_len - 2), 3):
                        stop = seq[j : j + 3]
                        if stop in STOP_CODONS:
                            orf_len = j - i
                            if min_len <= orf_len <= max_len:
                                orfs.append((i, j))
                            break
                i += 3

        return orfs

    @staticmethod
    def _overlaps_existing(start: int, end: int, occupied: list[tuple[int, int]]) -> bool:
        """Check if a region overlaps any existing feature."""
        for occ_start, occ_end in occupied:
            if start <= occ_end and end >= occ_start:
                return True
        return False
