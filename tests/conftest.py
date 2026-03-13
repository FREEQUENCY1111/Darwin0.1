"""Test fixtures for the Darwin ecosphere."""

from pathlib import Path

import pytest

from darwin.models import (
    Contig,
    Feature,
    FeatureType,
    Genome,
    Strand,
)
from darwin.rocks.models import (
    AnnotationConfig,
    Contig as RocksContig,
    Feature as RocksFeature,
    FeatureType as RocksFeatureType,
    Genome as RocksGenome,
    Strand as RocksStrand,
)
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Stream


@pytest.fixture
def sample_genome() -> Genome:
    """A minimal prokaryotic genome for testing (darwin.models API)."""
    seq = "ATGCGTACGATCGATCG" * 100  # 1700 bp
    return Genome(
        name="test_genome",
        contigs=[
            Contig(
                id="contig_1",
                sequence=seq,
                description="Test contig",
            ),
            Contig(
                id="contig_2",
                sequence=seq,
                description="Test contig 2",
            ),
        ],
    )


@pytest.fixture
def annotated_genome() -> Genome:
    """A genome with some features already attached (darwin.models API)."""
    seq = "ATGCGTACGATCGATCG" * 200  # 3400 bp
    genome = Genome(
        name="test_genome",
        contigs=[
            Contig(
                id="contig_1",
                sequence=seq,
                description="Annotated test contig",
                features=[
                    Feature(
                        seq_id="contig_1",
                        feature_type=FeatureType.CDS,
                        start=1,
                        end=900,
                        strand=Strand.FORWARD,
                        attributes={"locus_tag": "TEST_00001", "product": "test protein"},
                    ),
                    Feature(
                        seq_id="contig_1",
                        feature_type=FeatureType.TRNA,
                        start=1900,
                        end=1975,
                        strand=Strand.FORWARD,
                        attributes={"locus_tag": "TEST_t0001", "product": "tRNA-Ala"},
                    ),
                ],
            ),
            Contig(
                id="contig_2",
                sequence=seq,
                description="Test contig 2",
            ),
        ],
    )
    return genome


@pytest.fixture
def tmp_fasta(sample_genome: Genome, tmp_path: Path) -> Path:
    """Write a temporary FASTA file."""
    fasta_path = tmp_path / "test_genome.fasta"
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


# === Fixtures for darwin.rocks module tests ===
# These provide the rocks.models API instead of darwin.models API


@pytest.fixture
def rocks_sample_genome() -> RocksGenome:
    """A minimal prokaryotic genome for rocks testing (rocks.models API)."""
    seq = "ATGCGTACGATCGATCG" * 100  # 1700 bp
    return RocksGenome(
        name="test_genome",
        contigs=[
            RocksContig(
                id="contig_1",
                sequence=seq,
                description="Test contig",
            ),
        ],
    )


@pytest.fixture
def rocks_annotated_genome() -> RocksGenome:
    """A genome with features (rocks.models API)."""
    seq = "ATGCGTACGATCGATCG" * 200  # 3400 bp
    genome = RocksGenome(
        name="test_genome",
        contigs=[
            RocksContig(
                id="contig_1",
                sequence=seq,
                description="Annotated test contig",
                features=[
                    RocksFeature(
                        type=RocksFeatureType.CDS,
                        start=1,
                        end=900,
                        strand=RocksStrand.FORWARD,
                        contig_id="contig_1",
                        locus_tag="TEST_00001",
                        product="test protein",
                        translation="MRVDRS" * 50,
                    ),
                    RocksFeature(
                        type=RocksFeatureType.CDS,
                        start=950,
                        end=1800,
                        strand=RocksStrand.FORWARD,
                        contig_id="contig_1",
                        locus_tag="TEST_00002",
                        product="DNA gyrase subunit A",
                        translation="MKTLP" * 56,
                    ),
                ],
            )
        ],
    )
    return genome


@pytest.fixture
def rocks_tmp_fasta(rocks_sample_genome: RocksGenome, tmp_path: Path) -> Path:
    """Write a temporary FASTA file for rocks testing."""
    fasta_path = tmp_path / "test_genome.fasta"
    with open(fasta_path, "w") as fh:
        for contig in rocks_sample_genome.contigs:
            fh.write(f">{contig.id} {contig.description}\n")
            fh.write(contig.sequence + "\n")
    return fasta_path
