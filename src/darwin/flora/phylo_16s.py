"""
PhyloIdentifier — 16S-based taxonomic identification organism.

Feeds on: rrna.detected (rRNA features including 16S)
Produces: taxonomy.inferred (taxonomic classification)

Like a naturalist identifying species by their song —
uses 16S rRNA diagnostic probe sequences to identify the
organism's taxonomic group.

Algorithm:
  1. Extract 16S rRNA gene sequence (strand-aware)
  2. Domain classification: Bacteria vs Archaea
     (EUB338 target present → Bacteria; absent → Archaea)
  3. Phylum/class ID using known FISH probe target sites
     with fuzzy matching (up to 2 mismatches)
"""

from __future__ import annotations

import logging

from darwin.flora.base import Organism
from darwin.rocks.models import FeatureType, Genome, Strand
from darwin.soil.nutrients import NutrientStore
from darwin.soil.taxonomy_signatures import (
    ARCHAEAL_MARKERS,
    ARCHAEAL_PROBES,
    BACTERIAL_MARKERS,
    TAXONOMY_PROBES,
)
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.flora.phylo_16s")


# ── helpers ─────────────────────────────────────────────


def _reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence."""
    comp = str.maketrans("ATCGatcgNn", "TAGCtagcNn")
    return seq.translate(comp)[::-1]


def _extract_16s_sequences(genome: Genome) -> list[str]:
    """
    Extract all 16S rRNA sequences from the genome.

    Returns the *sense-strand* sequence (matching the rRNA)
    regardless of which DNA strand the gene is encoded on.
    """
    sequences = []
    rrnas = genome.features_of_type(FeatureType.RRNA)
    for f in rrnas:
        if "16S" in f.product:
            for contig in genome.contigs:
                if contig.id == f.contig_id:
                    seq = contig.sequence[f.start - 1 : f.end].upper()
                    # If gene is on reverse strand, flip so we always
                    # work with the rRNA-sense sequence
                    if f.strand == Strand.REVERSE:
                        seq = _reverse_complement(seq)
                    if len(seq) >= 500:
                        sequences.append(seq)
                    break
    return sequences


def _fuzzy_contains(haystack: str, needle: str, max_mismatches: int = 2) -> bool:
    """Check if *needle* appears in *haystack* with ≤ max_mismatches."""
    n = len(needle)
    if n > len(haystack):
        return False
    # Quick exact check first (common case)
    if needle in haystack:
        return True
    # Sliding-window fuzzy match
    for i in range(len(haystack) - n + 1):
        mismatches = 0
        for j in range(n):
            if haystack[i + j] != needle[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        else:
            return True
    return False


def _classify_domain(seq: str) -> str:
    """Classify as 'Bacteria', 'Archaea', or 'Unknown'."""
    # Check for bacterial markers (EUB338 target site)
    bact_hits = sum(
        1 for marker in BACTERIAL_MARKERS if _fuzzy_contains(seq, marker, 2)
    )
    # Check for archaeal markers
    arch_hits = sum(
        1 for marker in ARCHAEAL_MARKERS if _fuzzy_contains(seq, marker, 2)
    )

    if bact_hits > 0 and arch_hits == 0:
        return "Bacteria"
    if arch_hits > 0 and bact_hits == 0:
        return "Archaea"
    if arch_hits > bact_hits:
        return "Archaea"
    if bact_hits > arch_hits:
        return "Bacteria"

    # Fallback: sequence length heuristic
    # Bacterial 16S ≈ 1500-1550bp, Archaeal ≈ 1470-1500bp
    if len(seq) >= 1520:
        return "Bacteria"
    if len(seq) <= 1490:
        return "Archaea"
    return "Unknown"


def _classify_bacterial_group(seq: str) -> list:
    """Score each bacterial phylum/class against the 16S sequence."""
    results = []
    for group_name, info in TAXONOMY_PROBES.items():
        probes = info["probes"]
        max_mm = info.get("max_mismatches", 2)
        hits = sum(1 for p in probes if _fuzzy_contains(seq, p, max_mm))
        score = hits / len(probes) if probes else 0.0
        # Weight by number of probes that hit — more probes hitting
        # gives higher confidence than a single probe
        raw_hits = hits
        results.append((group_name, score, info["description"], raw_hits))

    # Sort by score first, then by raw hit count to break ties
    results.sort(key=lambda x: (x[1], x[3]), reverse=True)
    # Return in original (name, score, desc) format
    return [(name, score, desc) for name, score, desc, _ in results]


def _classify_archaeal_group(seq: str) -> list:
    """Score each archaeal phylum against the 16S sequence."""
    results = []
    for group_name, info in ARCHAEAL_PROBES.items():
        probes = info["probes"]
        max_mm = info.get("max_mismatches", 2)
        hits = sum(1 for p in probes if _fuzzy_contains(seq, p, max_mm))
        score = hits / len(probes) if probes else 0.0
        results.append((group_name, score, info["description"]))

    results.sort(key=lambda x: x[1], reverse=True)
    return results


def _score_taxonomy(seq: str) -> tuple[str, float, str, list]:
    """
    Hierarchical 16S taxonomy classification.

    Returns: (best_group, confidence, description, top_matches)
    """
    domain = _classify_domain(seq)

    if domain == "Bacteria":
        rankings = _classify_bacterial_group(seq)
        # Only accept if the top match actually had a probe hit
        if rankings and rankings[0][1] > 0:
            top_name, top_score, top_desc = rankings[0]
            top_matches = [
                {"group": name, "score": round(score, 3), "description": desc}
                for name, score, desc in rankings[:3]
                if score > 0
            ]
            return top_name, top_score, top_desc, top_matches
        else:
            # Domain known, but no phylum-level match
            return "Bacteria", 0.5, "Domain Bacteria (phylum uncertain)", []

    elif domain == "Archaea":
        rankings = _classify_archaeal_group(seq)
        if rankings and rankings[0][1] > 0:
            top_name, top_score, top_desc = rankings[0]
            top_matches = [
                {"group": name, "score": round(score, 3), "description": desc}
                for name, score, desc in rankings[:3]
                if score > 0
            ]
            return top_name, top_score, top_desc, top_matches
        else:
            return "Archaea", 0.5, "Domain Archaea (phylum uncertain)", []

    else:
        return "", 0.0, "", []


class PhyloIdentifier(Organism):
    """16S phylogenetic identifier — the naturalist."""

    name = "phylo_16s"
    feeds_on_nutrients = [NutrientType.RRNA_DETECTED]
    produces_nutrients = [NutrientType.TAXONOMY_INFERRED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return True  # Pure Python

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """Identify organism taxonomy from 16S rRNA."""
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})

        self.logger.info("🔬 Identifying taxonomy from 16S rRNA...")

        sequences = _extract_16s_sequences(genome)

        if not sequences:
            self.logger.info("No 16S rRNA sequences found — skipping taxonomy")
            return Nutrient(
                type=NutrientType.TAXONOMY_INFERRED,
                data={
                    "genome": genome,
                    "taxonomy": "",
                    "confidence": 0.0,
                    "top_matches": [],
                    "config": config,
                },
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        # Use the longest 16S sequence for identification
        best_seq = max(sequences, key=len)

        # Hierarchical classification
        top_name, confidence, description, top_matches = _score_taxonomy(best_seq)

        if top_name:
            genome.taxonomy = top_name
            if not genome.organism:
                genome.organism = f"unclassified {top_name}"

            self.logger.info(
                f"🔬 Taxonomy: {top_name} — {description} "
                f"(confidence: {confidence:.1%})"
            )
        else:
            self.logger.info("🔬 Could not determine taxonomy from 16S")

        return Nutrient(
            type=NutrientType.TAXONOMY_INFERRED,
            data={
                "genome": genome,
                "taxonomy": top_name,
                "confidence": round(confidence, 3),
                "top_matches": top_matches,
                "num_16s_sequences": len(sequences),
                "best_16s_length": len(best_seq),
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )
