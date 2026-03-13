"""Tests for the Rocks layer — data models and FASTA I/O."""

from pathlib import Path

import pytest

from darwin.rocks.fasta import parse_fasta, write_fasta, write_proteins
from darwin.rocks.models import Feature, FeatureType, Strand


class TestGenome:
    def test_total_length(self, rocks_sample_genome):
        assert rocks_sample_genome.total_length == 1700

    def test_gc_content(self, rocks_sample_genome):
        gc = rocks_sample_genome.gc_content
        assert 0 < gc < 100

    def test_num_contigs(self, rocks_sample_genome):
        assert rocks_sample_genome.num_contigs == 1

    def test_summary(self, rocks_sample_genome):
        s = rocks_sample_genome.summary()
        assert s["name"] == "test_genome"
        assert s["num_contigs"] == 1
        assert s["total_length_bp"] == 1700


class TestFeature:
    def test_length(self):
        f = Feature(type=FeatureType.CDS, start=100, end=400)
        assert f.length == 300

    def test_location_str_forward(self):
        f = Feature(type=FeatureType.CDS, start=1, end=900, strand=Strand.FORWARD)
        assert f.location_str == "1..900"

    def test_location_str_reverse(self):
        f = Feature(type=FeatureType.CDS, start=1, end=900, strand=Strand.REVERSE)
        assert "complement" in f.location_str


class TestFasta:
    def test_parse_fasta(self, rocks_tmp_fasta):
        genome = parse_fasta(rocks_tmp_fasta)
        assert genome.num_contigs == 1
        assert genome.total_length > 0

    def test_parse_missing_file(self):
        with pytest.raises(FileNotFoundError):
            parse_fasta(Path("/nonexistent.fasta"))

    def test_write_fasta(self, rocks_sample_genome, tmp_path):
        out = tmp_path / "out.fasta"
        result = write_fasta(rocks_sample_genome, out)
        assert result.exists()
        content = result.read_text()
        assert ">contig_1" in content

    def test_write_proteins(self, rocks_annotated_genome, tmp_path):
        out = tmp_path / "proteins.faa"
        result = write_proteins(rocks_annotated_genome, out)
        assert result.exists()
        content = result.read_text()
        assert ">TEST_00001" in content
