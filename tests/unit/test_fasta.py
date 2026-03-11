"""Tests for FASTA I/O utilities."""

from pathlib import Path

from darwin.utils.fasta import parse_fasta, write_fasta


class TestParseFasta:
    def test_parse_basic(self, tmp_fasta: Path):
        genome = parse_fasta(tmp_fasta, min_len=100)
        assert genome.num_contigs == 2
        assert genome.name == "test_genome"

    def test_parse_filters_short(self, tmp_path: Path):
        fasta = tmp_path / "short.fasta"
        fasta.write_text(">short\nATG\n>long\n" + "ATGC" * 100 + "\n")
        genome = parse_fasta(fasta, min_len=200)
        assert genome.num_contigs == 1
        assert genome.contigs[0].id == "long"

    def test_parse_empty(self, tmp_path: Path):
        fasta = tmp_path / "empty.fasta"
        fasta.write_text("")
        genome = parse_fasta(fasta)
        assert genome.num_contigs == 0


class TestWriteFasta:
    def test_roundtrip(self, tmp_path: Path, sample_genome):
        out = tmp_path / "out.fasta"
        write_fasta(sample_genome.contigs, out)
        assert out.exists()
        # Re-parse and verify
        genome2 = parse_fasta(out, min_len=0)
        assert genome2.num_contigs == 2
        assert genome2.contigs[0].id == "contig_1"
