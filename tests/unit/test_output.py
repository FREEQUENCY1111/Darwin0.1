"""Tests for output format writers."""

import json
from pathlib import Path

from darwin.models import Genome
from darwin.output.gff import write_gff3
from darwin.output.genbank import write_genbank
from darwin.output.json_out import genome_to_dict, write_json


class TestGFF3:
    def test_write_empty_genome(self, tmp_path: Path, sample_genome: Genome):
        out = tmp_path / "test.gff3"
        result = write_gff3(sample_genome, out)
        assert result.exists()
        content = out.read_text()
        assert "##gff-version 3" in content

    def test_write_with_features(self, tmp_path: Path, annotated_genome: Genome):
        out = tmp_path / "test.gff3"
        write_gff3(annotated_genome, out)
        content = out.read_text()
        assert "CDS" in content
        assert "tRNA" in content
        assert "TEST_00001" in content


class TestGenBank:
    def test_write_empty_genome(self, tmp_path: Path, sample_genome: Genome):
        out = tmp_path / "test.gbk"
        result = write_genbank(sample_genome, out)
        assert result.exists()

    def test_write_with_features(self, tmp_path: Path, annotated_genome: Genome):
        out = tmp_path / "test.gbk"
        write_genbank(annotated_genome, out)
        content = out.read_text()
        assert "LOCUS" in content
        assert "test protein" in content


class TestJSON:
    def test_write(self, tmp_path: Path, annotated_genome: Genome):
        out = tmp_path / "test.json"
        write_json(annotated_genome, out)
        data = json.loads(out.read_text())
        assert data["genome"]["name"] == "test_genome"
        assert len(data["contigs"]) == 2

    def test_genome_to_dict(self, annotated_genome: Genome):
        d = genome_to_dict(annotated_genome)
        assert d["genome"]["cds_count"] == 1
        assert d["genome"]["trna_count"] == 1
