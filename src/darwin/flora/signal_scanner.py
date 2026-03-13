"""
SignalScanner — Signal peptide detection organism.

Feeds on: genes.called (predicted CDS with translations)
Produces: signal_peptides.found

Like a customs inspector checking each protein for an
export tag — signal peptides are N-terminal sequences
that direct proteins to the cell membrane or secretion.

Algorithm: for each CDS translation, check first 70 AA:
1. Positively charged n-region (first 5 AA)
2. Hydrophobic h-region (Kyte-Doolittle window ≥ 1.6)
3. Cleavage site [A/G]-X-[A/S] at positions 15-30
"""

from __future__ import annotations

import logging

from darwin.flora.base import Organism
from darwin.rocks.models import Feature, FeatureType, Genome, Strand
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.flora.signal_scanner")

# Kyte-Doolittle hydropathy scale
KD_SCALE: dict[str, float] = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "D": -3.5,
    "E": -3.5, "N": -3.5, "Q": -3.5, "K": -3.9, "R": -4.5,
}

# Positive charge residues
POSITIVE_AA = {"K", "R", "H"}

# Cleavage site: small/neutral residues at -1 and -3 (von Heijne rule)
CLEAVAGE_M3 = {"A", "G", "S", "T", "V", "L", "I"}  # -3 position
CLEAVAGE_M1 = {"A", "G", "S", "T"}  # -1 position


def _hydropathy_window(seq: str, window: int = 8) -> list[float]:
    """Compute Kyte-Doolittle hydropathy over sliding window."""
    scores = []
    for i in range(len(seq) - window + 1):
        window_seq = seq[i : i + window]
        avg = sum(KD_SCALE.get(aa, 0.0) for aa in window_seq) / window
        scores.append(avg)
    return scores


def _detect_signal_peptide(translation: str) -> dict | None:
    """
    Check if a protein has a signal peptide.

    Returns dict with cleavage_site position and score, or None.
    """
    if len(translation) < 30:
        return None

    # Work with first 70 AA (signal peptides are always N-terminal)
    sp_region = translation[:70].upper()

    # 1. Check n-region: first 5 residues should have positive charge
    n_region = sp_region[:5]
    positive_count = sum(1 for aa in n_region if aa in POSITIVE_AA)
    if positive_count < 1:
        return None

    # 2. Check h-region: hydrophobic core in positions 2-20
    h_region = sp_region[1:20]
    hydro_scores = _hydropathy_window(h_region, window=8)
    if not hydro_scores:
        return None

    max_hydro = max(hydro_scores)
    if max_hydro < 1.6:
        return None

    # Find h-region boundaries
    h_start = 1
    h_end = 1
    for i, score in enumerate(hydro_scores):
        if score >= 1.0:
            h_start = i + 1
            break
    for i in range(len(hydro_scores) - 1, -1, -1):
        if hydro_scores[i] >= 1.0:
            h_end = i + 1 + 8  # end of window
            break

    # 3. Check cleavage site: [A/G/S/T/V/L/I]-X-[A/G/S/T] at positions 15-35
    best_site = None
    best_score = 0.0

    search_start = max(15, h_end - 2)
    search_end = min(35, len(sp_region))

    for pos in range(search_start, search_end):
        if pos < 2:
            continue

        # Von Heijne's (-3, -1) rule
        m3 = sp_region[pos - 3] if pos >= 3 else ""
        m1 = sp_region[pos - 1]

        if m3 in CLEAVAGE_M3 and m1 in CLEAVAGE_M1:
            # Score: higher is better
            score = 0.0
            score += 1.0 if m3 in {"A", "G"} else 0.5
            score += 1.0 if m1 in {"A", "G"} else 0.5
            score += max_hydro / 4.5  # normalize hydropathy contribution

            if score > best_score:
                best_score = score
                best_site = pos

    if best_site is None:
        return None

    return {
        "cleavage_site": best_site,
        "score": round(best_score, 2),
        "n_positive": positive_count,
        "h_max": round(max_hydro, 2),
        "h_region": f"{h_start}-{h_end}",
    }


class SignalScanner(Organism):
    """Signal peptide detector — the customs inspector."""

    name = "signal_scanner"
    feeds_on_nutrients = [NutrientType.GENES_CALLED]
    produces_nutrients = [NutrientType.SIGNAL_PEPTIDES_FOUND]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return True  # Pure Python, no external tools needed

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """Scan all CDS translations for signal peptides."""
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})

        self.logger.info("📮 Scanning for signal peptides...")

        cds_features = genome.features_of_type(FeatureType.CDS)
        sp_count = 0

        for f in cds_features:
            if not f.translation:
                continue

            result = _detect_signal_peptide(f.translation)
            if result:
                sp_count += 1
                # Annotate the CDS with signal peptide info
                sp_note = (
                    f"signal_peptide:1..{result['cleavage_site']};"
                    f"score={result['score']};"
                    f"h_max={result['h_max']}"
                )
                if f.note:
                    f.note += f";{sp_note}"
                else:
                    f.note = sp_note

                # Add a signal_peptide feature to the contig
                for contig in genome.contigs:
                    if contig.id == f.contig_id:
                        # Calculate signal peptide genomic coordinates
                        sp_aa_len = result["cleavage_site"]
                        sp_nt_len = sp_aa_len * 3

                        if f.strand == Strand.FORWARD:
                            sp_start = f.start
                            sp_end = f.start + sp_nt_len - 1
                        else:
                            sp_end = f.end
                            sp_start = f.end - sp_nt_len + 1

                        sp_feature = Feature(
                            type=FeatureType.SIGNAL_PEPTIDE,
                            start=sp_start,
                            end=sp_end,
                            strand=f.strand,
                            score=result["score"],
                            contig_id=contig.id,
                            locus_tag="",
                            product=f"signal peptide ({f.locus_tag})",
                            inference="ab initio prediction:Darwin:SignalScanner",
                            note=f"parent={f.locus_tag};cleavage_site={result['cleavage_site']}",
                        )
                        contig.features.append(sp_feature)
                        break

        self.logger.info(
            f"📮 Found {sp_count} signal peptides in "
            f"{len(cds_features)} proteins"
        )

        return Nutrient(
            type=NutrientType.SIGNAL_PEPTIDES_FOUND,
            data={
                "genome": genome,
                "signal_peptide_count": sp_count,
                "total_cds": len(cds_features),
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )
