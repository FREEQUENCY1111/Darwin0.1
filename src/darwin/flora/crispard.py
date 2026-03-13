"""
CRISPARd — CRISPR array detection organism.

Feeds on: genome.loaded (raw sequence data)
Produces: crispr.detected (CRISPR arrays with repeats and spacers)

Like a sentinel organism that scans the environment for signs
of ancient viral warfare — CRISPR arrays are molecular fossils
of past phage infections.

Algorithm: repeat-spacer pattern matching — find direct repeats
23-47bp with spacers 26-72bp between them, requiring ≥3 units.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from darwin.flora.base import Organism
from darwin.rocks.models import Feature, FeatureType, Genome, Strand
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.flora.crispard")

# CRISPR repeat/spacer constraints
MIN_REPEAT_LEN = 23
MAX_REPEAT_LEN = 47
MIN_SPACER_LEN = 26
MAX_SPACER_LEN = 72
MIN_UNITS = 3  # minimum repeat-spacer units to call a CRISPR array
MAX_MISMATCH_RATIO = 0.15  # allow 15% mismatches between repeat copies


@dataclass
class CRISPRArray:
    """A detected CRISPR array."""

    contig_id: str
    start: int
    end: int
    repeat_seq: str
    repeat_len: int
    spacer_count: int
    spacers: list[str]
    repeat_copies: int
    consensus_score: float  # how conserved the repeats are (0-1)


def _hamming_distance(s1: str, s2: str) -> int:
    """Count mismatches between two equal-length strings."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2, strict=False))


def _find_crispr_arrays(sequence: str, contig_id: str) -> list[CRISPRArray]:
    """
    Find CRISPR arrays in a sequence using repeat-spacer pattern matching.

    Strategy:
    1. For each possible repeat length (23-47bp), scan for exact k-mer repeats
    2. Check if repeats are spaced at regular intervals with spacers in range
    3. Extend arrays by finding additional repeat copies with mismatches
    4. Keep arrays with ≥3 repeat-spacer units
    """
    seq = sequence.upper()
    seq_len = len(seq)
    arrays: list[CRISPRArray] = []

    if seq_len < (MIN_REPEAT_LEN + MIN_SPACER_LEN) * MIN_UNITS:
        return arrays

    # Track regions already assigned to arrays
    used_regions: list[tuple[int, int]] = []

    # Scan across repeat lengths
    for repeat_len in range(MIN_REPEAT_LEN, min(MAX_REPEAT_LEN + 1, seq_len // 3)):
        # Build k-mer index for this repeat length
        kmer_positions: dict[str, list[int]] = {}
        for i in range(seq_len - repeat_len + 1):
            kmer = seq[i : i + repeat_len]
            if "N" in kmer:
                continue
            if kmer not in kmer_positions:
                kmer_positions[kmer] = []
            kmer_positions[kmer].append(i)

        # Find k-mers with multiple occurrences
        for kmer, positions in kmer_positions.items():
            if len(positions) < MIN_UNITS:
                continue

            # Check if positions form a regular array pattern
            positions.sort()
            array = _try_build_array(seq, positions, kmer, repeat_len, contig_id)
            if array and not _overlaps_used(array.start, array.end, used_regions):
                arrays.append(array)
                used_regions.append((array.start, array.end))

    # Sort by position
    arrays.sort(key=lambda a: a.start)
    return arrays


def _try_build_array(
    seq: str,
    positions: list[int],
    repeat: str,
    repeat_len: int,
    contig_id: str,
) -> CRISPRArray | None:
    """Try to build a CRISPR array from repeat positions."""
    # Check consecutive positions for valid spacer lengths
    units = [positions[0]]

    for i in range(1, len(positions)):
        gap = positions[i] - positions[i - 1] - repeat_len
        if MIN_SPACER_LEN <= gap <= MAX_SPACER_LEN:
            if not units or positions[i - 1] == units[-1]:
                if not units:
                    units.append(positions[i - 1])
                units.append(positions[i])
        else:
            # Break in pattern — check if we have enough units
            if len(units) >= MIN_UNITS:
                break
            units = [positions[i]]

    if len(units) < MIN_UNITS:
        return None

    # Extract spacers
    spacers = []
    for i in range(len(units) - 1):
        spacer_start = units[i] + repeat_len
        spacer_end = units[i + 1]
        spacers.append(seq[spacer_start:spacer_end])

    # Compute repeat conservation score
    repeats = [seq[pos : pos + repeat_len] for pos in units]
    if len(repeats) > 1:
        total_mismatches = 0
        comparisons = 0
        for i in range(1, len(repeats)):
            total_mismatches += _hamming_distance(repeats[0], repeats[i])
            comparisons += 1
        avg_mismatch_rate = total_mismatches / (comparisons * repeat_len) if comparisons else 0
        consensus_score = 1.0 - avg_mismatch_rate
    else:
        consensus_score = 1.0

    # Only accept if repeats are well conserved
    if consensus_score < (1.0 - MAX_MISMATCH_RATIO):
        return None

    array_start = units[0] + 1  # 1-based
    array_end = units[-1] + repeat_len

    return CRISPRArray(
        contig_id=contig_id,
        start=array_start,
        end=array_end,
        repeat_seq=repeat,
        repeat_len=repeat_len,
        spacer_count=len(spacers),
        spacers=spacers,
        repeat_copies=len(units),
        consensus_score=consensus_score,
    )


def _overlaps_used(start: int, end: int, used: list[tuple[int, int]]) -> bool:
    for us, ue in used:
        if start <= ue and end >= us:
            return True
    return False


class CRISPARd(Organism):
    """CRISPR array detector — the sentinel."""

    name = "crispard"
    feeds_on_nutrients = [NutrientType.GENOME_LOADED]
    produces_nutrients = [NutrientType.CRISPR_DETECTED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return True  # Pure Python, no external tools needed

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """Scan genome for CRISPR arrays."""
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})

        self.logger.info("🛡️ Scanning for CRISPR arrays...")

        total_arrays = 0
        total_spacers = 0

        for contig in genome.contigs:
            arrays = _find_crispr_arrays(contig.sequence, contig.id)

            for array in arrays:
                # Add CRISPR feature to contig
                feature = Feature(
                    type=FeatureType.CRISPR,
                    start=array.start,
                    end=array.end,
                    strand=Strand.UNKNOWN,
                    score=array.consensus_score,
                    contig_id=contig.id,
                    locus_tag="",
                    product=f"CRISPR array ({array.repeat_copies} repeats, {array.spacer_count} spacers)",
                    inference="ab initio prediction:Darwin:CRISPARd",
                    note=(
                        f"repeat_len={array.repeat_len};"
                        f"spacer_count={array.spacer_count};"
                        f"consensus={array.consensus_score:.2f};"
                        f"repeat_seq={array.repeat_seq[:30]}"
                    ),
                )
                contig.features.append(feature)
                total_arrays += 1
                total_spacers += array.spacer_count

        self.logger.info(
            f"🛡️ Found {total_arrays} CRISPR array(s) "
            f"with {total_spacers} total spacers"
        )

        return Nutrient(
            type=NutrientType.CRISPR_DETECTED,
            data={
                "genome": genome,
                "array_count": total_arrays,
                "total_spacers": total_spacers,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )
