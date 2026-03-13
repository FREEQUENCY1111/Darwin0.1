"""
ANI — Average Nucleotide Identity estimator.

Approximates ANI between two genomes using k-mer Jaccard
similarity. Like comparing rock strata — similar mineral
composition suggests shared geological history.

Algorithm:
  1. Compute k-mer sets (k=16) for both genomes
  2. Calculate Jaccard similarity = |A ∩ B| / |A ∪ B|
  3. Convert to approximate ANI using Mash-style formula:
     ANI ≈ 1 + (1/k) * ln(2J / (1+J))

Species boundary: ANI ≥ 95% (Goris et al., 2007)
"""

from __future__ import annotations

import logging
from pathlib import Path

from darwin.rocks.fasta import parse_fasta

logger = logging.getLogger("darwin.utils.ani")

DEFAULT_K = 16
SKETCH_SIZE = 10000  # max k-mers to keep per genome (MinHash-style)


def compute_kmer_set(sequence: str, k: int = DEFAULT_K) -> set[str]:
    """Extract all k-mers from a sequence, both strands."""
    kmers: set[str] = set()
    seq = sequence.upper()
    rc_map = str.maketrans("ACGT", "TGCA")

    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" not in kmer:
            # Canonical k-mer (min of forward and reverse complement)
            rc = kmer.translate(rc_map)[::-1]
            kmers.add(min(kmer, rc))

    return kmers


def sketch_genome(sequence: str, k: int = DEFAULT_K, sketch_size: int = SKETCH_SIZE) -> set[str]:
    """
    Create a MinHash-style sketch of a genome.

    Keeps only the `sketch_size` lexicographically smallest k-mers
    for memory-efficient comparison of large genomes.
    """
    all_kmers = compute_kmer_set(sequence, k)

    if len(all_kmers) <= sketch_size:
        return all_kmers

    # Keep smallest k-mers (deterministic, reproducible)
    sorted_kmers = sorted(all_kmers)
    return set(sorted_kmers[:sketch_size])


def jaccard_similarity(set_a: set[str], set_b: set[str]) -> float:
    """Compute Jaccard similarity between two sets."""
    if not set_a or not set_b:
        return 0.0
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union > 0 else 0.0


def jaccard_to_ani(jaccard: float, k: int = DEFAULT_K) -> float:
    """
    Convert Jaccard similarity to approximate ANI.

    Uses the Mash formula (Ondov et al., 2016):
    ANI ≈ 1 + (1/k) * ln(2J / (1+J))
    """
    import math

    if jaccard <= 0:
        return 0.0
    if jaccard >= 1.0:
        return 100.0

    mash_distance = -1.0 / k * math.log(2.0 * jaccard / (1.0 + jaccard))
    ani = (1.0 - mash_distance) * 100.0

    return max(0.0, min(100.0, round(ani, 2)))


def compare_genomes(
    genome1_path: Path,
    genome2_path: Path,
    k: int = DEFAULT_K,
    sketch_size: int = SKETCH_SIZE,
) -> dict:
    """
    Compare two genomes and estimate ANI.

    Returns a dict with Jaccard similarity, approximate ANI,
    shared k-mer count, and species-level determination.
    """
    logger.info(f"🔬 Comparing {genome1_path.name} vs {genome2_path.name}...")

    # Parse genomes
    g1 = parse_fasta(genome1_path)
    g2 = parse_fasta(genome2_path)

    # Concatenate all contigs
    seq1 = "".join(c.sequence for c in g1.contigs)
    seq2 = "".join(c.sequence for c in g2.contigs)

    logger.info(f"  Genome 1: {g1.name} ({len(seq1):,} bp, {g1.num_contigs} contigs)")
    logger.info(f"  Genome 2: {g2.name} ({len(seq2):,} bp, {g2.num_contigs} contigs)")

    # Sketch both genomes
    sketch1 = sketch_genome(seq1, k, sketch_size)
    sketch2 = sketch_genome(seq2, k, sketch_size)

    # Compute similarity
    jaccard = jaccard_similarity(sketch1, sketch2)
    ani = jaccard_to_ani(jaccard, k)
    shared = len(sketch1 & sketch2)

    same_species = ani >= 95.0

    result = {
        "genome1": {
            "name": g1.name,
            "length_bp": len(seq1),
            "contigs": g1.num_contigs,
            "gc_content": g1.gc_content,
            "kmers": len(sketch1),
        },
        "genome2": {
            "name": g2.name,
            "length_bp": len(seq2),
            "contigs": g2.num_contigs,
            "gc_content": g2.gc_content,
            "kmers": len(sketch2),
        },
        "comparison": {
            "k": k,
            "sketch_size": sketch_size,
            "jaccard_similarity": round(jaccard, 6),
            "approximate_ani": ani,
            "shared_kmers": shared,
            "same_species": same_species,
            "threshold": "ANI >= 95% (Goris et al., 2007)",
        },
    }

    logger.info(f"  Jaccard: {jaccard:.4f}, ANI: {ani:.2f}%")
    if same_species:
        logger.info("  ✅ Same species (ANI ≥ 95%)")
    else:
        logger.info("  ❌ Different species (ANI < 95%)")

    return result
