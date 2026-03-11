"""Shared test fixtures for Darwin."""

from __future__ import annotations

from pathlib import Path

import pytest

from darwin.models import (
    AnnotationConfig,
    Contig,
    Feature,
    FeatureType,
    Genome,
    Strand,
)

# Small test genome — 2 contigs of realistic prokaryotic sequence
# This is a truncated snippet, enough for unit tests
TEST_SEQ_1 = (
    "ATGAAAGCGATTATTGGTCTGGGTGCTTATAAT"
    "CTGACCGCAGAAGATCTGACCGAAATCAACGAC"
    "TGGAAAGCGGATAACACCCCGCTGACCCTGCGT"
    "GATCTGTCAGAAGCGCGTCGTCTGTGGGCACTG"
    "GCAACCGATGATATTCTGAAAGGCCTGATCGGT"
    "TTTCCGGCAAGCGAGCATCTGATCGTGTAA"
) * 20  # ~3.6kb

TEST_SEQ_2 = (
    "GTGCGTAAAGCGATCGTTCCGGGTTTTACCGAT"
    "ACCCTGCCAGCAGATCTGACCGAAGTGGAAGCG"
    "TGGCGTGCGGATAACGGCCCGCTGACCCTGGCG"
    "GATCTGTCAGAAGCGCGTCGTCTGTGGGCGCTG"
    "AAAGCGGATGATATTCTGAAAGGCCTGATCTAT"
    "TTTCCGGCAAGCGAGCATCTGATCATCTAA"
) * 15  # ~2.7kb


@pytest.fixture
def sample_genome() -> Genome:
    """A minimal test genome with 2 contigs."""
    return Genome(
        name="test_genome",
        contigs=[
            Contig(id="contig_1", sequence=TEST_SEQ_1, description="test contig 1"),
            Contig(id="contig_2", sequence=TEST_SEQ_2, description="test contig 2"),
        ],
    )


@pytest.fixture
def annotated_genome(sample_genome: Genome) -> Genome:
    """A genome with some pre-existing features for testing."""
    sample_genome.contigs[0].features = [
        Feature(
            seq_id="contig_1",
            feature_type=FeatureType.CDS,
            start=1,
            end=300,
            strand=Strand.FORWARD,
            attributes={"locus_tag": "TEST_00001", "product": "test protein"},
        ),
        Feature(
            seq_id="contig_1",
            feature_type=FeatureType.TRNA,
            start=400,
            end=475,
            strand=Strand.FORWARD,
            attributes={"locus_tag": "TEST_t0001", "product": "tRNA-Ala"},
        ),
    ]
    return sample_genome


@pytest.fixture
def tmp_fasta(tmp_path: Path, sample_genome: Genome) -> Path:
    """Write a temporary FASTA file from the sample genome."""
    fasta = tmp_path / "test_genome.fasta"
    with open(fasta, "w") as fh:
        for contig in sample_genome.contigs:
            fh.write(f">{contig.id} {contig.description}\n")
            seq = contig.sequence
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")
    return fasta


@pytest.fixture
def annotation_config(tmp_path: Path, tmp_fasta: Path) -> AnnotationConfig:
    """A test annotation config."""
    return AnnotationConfig(
        input_path=tmp_fasta,
        output_dir=tmp_path / "output",
        locus_tag_prefix="TEST",
    )
