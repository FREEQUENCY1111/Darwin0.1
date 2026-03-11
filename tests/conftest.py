"""Test fixtures for the Darwin ecosphere."""

import pytest
import tempfile
from pathlib import Path

from darwin.rocks.models import (
    Genome, Contig, Feature, FeatureType, Strand, AnnotationConfig,
)
from darwin.water.stream import Stream
from darwin.soil.nutrients import NutrientStore


@pytest.fixture
def sample_genome() -> Genome:
    """A minimal prokaryotic genome for testing."""
    seq = "ATGCGTACGATCGATCG" * 100  # 1700 bp
    return Genome(
        name="test_genome",
        contigs=[
            Contig(
                id="contig_1",
                sequence=seq,
                description="Test contig",
            )
        ],
    )


@pytest.fixture
def annotated_genome() -> Genome:
    """A genome with some features already attached."""
    seq = "ATGCGTACGATCGATCG" * 200  # 3400 bp
    genome = Genome(
        name="annotated_test",
        contigs=[
            Contig(
                id="contig_1",
                sequence=seq,
                description="Annotated test contig",
                features=[
                    Feature(
                        type=FeatureType.CDS,
                        start=1,
                        end=900,
                        strand=Strand.FORWARD,
                        contig_id="contig_1",
                        locus_tag="TEST_00001",
                        product="hypothetical protein",
                        translation="MRVDRS" * 50,
                    ),
                    Feature(
                        type=FeatureType.CDS,
                        start=950,
                        end=1800,
                        strand=Strand.FORWARD,
                        contig_id="contig_1",
                        locus_tag="TEST_00002",
                        product="DNA gyrase subunit A",
                        translation="MKTLP" * 56,
                    ),
                    Feature(
                        type=FeatureType.TRNA,
                        start=1900,
                        end=1975,
                        strand=Strand.FORWARD,
                        contig_id="contig_1",
                        locus_tag="TEST_t0001",
                        product="tRNA-Ala",
                    ),
                    Feature(
                        type=FeatureType.RRNA,
                        start=2000,
                        end=3540,
                        strand=Strand.FORWARD,
                        contig_id="contig_1",
                        locus_tag="TEST_r0001",
                        product="16S ribosomal RNA",
                    ),
                ],
            )
        ],
    )
    return genome


@pytest.fixture
def tmp_fasta(sample_genome: Genome, tmp_path: Path) -> Path:
    """Write a temporary FASTA file."""
    fasta_path = tmp_path / "test.fasta"
    with open(fasta_path, "w") as fh:
        for contig in sample_genome.contigs:
            fh.write(f">{contig.id} {contig.description}\n")
            fh.write(contig.sequence + "\n")
    return fasta_path


@pytest.fixture
def stream() -> Stream:
    """A fresh water stream."""
    return Stream()


@pytest.fixture
def soil() -> NutrientStore:
    """Barren soil (no tools installed in test env)."""
    return NutrientStore()


@pytest.fixture
def annotation_config(tmp_fasta: Path, tmp_path: Path) -> AnnotationConfig:
    """Basic annotation configuration."""
    return AnnotationConfig(
        input_file=tmp_fasta,
        output_dir=tmp_path / "output",
        locus_tag_prefix="TEST",
    )
