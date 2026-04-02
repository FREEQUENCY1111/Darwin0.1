"""Tests for the darwin domains feature — protein domain characterization."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from darwin.sunlight.domains import (
    DomainHit,
    ProteinDomains,
    format_all_maps,
    format_domain_map,
    format_tsv,
    read_protein_fasta,
)


# ── FASTA reader tests ──────────────────────────────────────────────────


class TestReadProteinFasta:
    def test_reads_single_protein(self, tmp_path: Path) -> None:
        fasta = tmp_path / "test.faa"
        fasta.write_text(">protein_1 some description\nMKTLPQRSTAA\nVKLMNOP\n")
        result = read_protein_fasta(fasta)
        assert result == {"protein_1": "MKTLPQRSTAAVKLMNOP"}

    def test_reads_multiple_proteins(self, tmp_path: Path) -> None:
        fasta = tmp_path / "test.faa"
        fasta.write_text(
            ">prot_A\nMKTLP\n"
            ">prot_B extra info\nVKLMN\nOPQRS\n"
            ">prot_C\nAAAA\n"
        )
        result = read_protein_fasta(fasta)
        assert len(result) == 3
        assert result["prot_A"] == "MKTLP"
        assert result["prot_B"] == "VKLMNOPQRS"
        assert result["prot_C"] == "AAAA"

    def test_empty_file(self, tmp_path: Path) -> None:
        fasta = tmp_path / "empty.faa"
        fasta.write_text("")
        result = read_protein_fasta(fasta)
        assert result == {}

    def test_takes_first_word_as_id(self, tmp_path: Path) -> None:
        fasta = tmp_path / "test.faa"
        fasta.write_text(">my_protein [Escherichia coli]\nMKTLP\n")
        result = read_protein_fasta(fasta)
        assert "my_protein" in result


# ── DomainHit model tests ───────────────────────────────────────────────


class TestDomainHit:
    def _make_hit(self, **kwargs) -> DomainHit:
        defaults = dict(
            protein_id="prot_1",
            protein_length=300,
            domain_name="PF00001",
            accession="PF00001.21",
            description="7tm_1",
            env_from=10,
            env_to=200,
            evalue=1e-30,
            score=100.5,
            db_name="Pfam-A",
        )
        defaults.update(kwargs)
        return DomainHit(**defaults)

    def test_domain_length(self) -> None:
        hit = self._make_hit(env_from=10, env_to=200)
        assert hit.domain_length == 191

    def test_coverage(self) -> None:
        hit = self._make_hit(env_from=1, env_to=150, protein_length=300)
        assert hit.coverage == pytest.approx(0.5, abs=0.01)

    def test_coverage_zero_length_protein(self) -> None:
        hit = self._make_hit(protein_length=0)
        assert hit.coverage == 0.0


# ── ProteinDomains model tests ──────────────────────────────────────────


class TestProteinDomains:
    def test_no_domains(self) -> None:
        pd = ProteinDomains(protein_id="prot_1", protein_length=300, hits=[])
        assert pd.num_domains == 0
        assert pd.total_coverage == 0.0

    def test_single_domain(self) -> None:
        hit = DomainHit(
            protein_id="prot_1", protein_length=300,
            domain_name="PF00001", accession="PF00001",
            description="test", env_from=1, env_to=150,
            evalue=1e-30, score=100.0, db_name="Pfam-A",
        )
        pd = ProteinDomains(protein_id="prot_1", protein_length=300, hits=[hit])
        assert pd.num_domains == 1
        assert pd.total_coverage == pytest.approx(0.5, abs=0.01)

    def test_two_non_overlapping_domains(self) -> None:
        hit1 = DomainHit(
            protein_id="prot_1", protein_length=400,
            domain_name="PF00001", accession="PF00001",
            description="", env_from=1, env_to=100,
            evalue=1e-30, score=100.0, db_name="Pfam-A",
        )
        hit2 = DomainHit(
            protein_id="prot_1", protein_length=400,
            domain_name="PF00002", accession="PF00002",
            description="", env_from=200, env_to=350,
            evalue=1e-20, score=80.0, db_name="Pfam-A",
        )
        pd = ProteinDomains(protein_id="prot_1", protein_length=400, hits=[hit1, hit2])
        assert pd.num_domains == 2
        # Coverage: (100 + 151) / 400 = 0.6275
        assert pd.total_coverage == pytest.approx(0.6275, abs=0.01)

    def test_overlapping_domains_merged(self) -> None:
        hit1 = DomainHit(
            protein_id="prot_1", protein_length=200,
            domain_name="PF00001", accession="PF00001",
            description="", env_from=1, env_to=120,
            evalue=1e-30, score=100.0, db_name="Pfam-A",
        )
        hit2 = DomainHit(
            protein_id="prot_1", protein_length=200,
            domain_name="PF00002", accession="PF00002",
            description="", env_from=100, env_to=200,
            evalue=1e-20, score=80.0, db_name="Pfam-A",
        )
        pd = ProteinDomains(protein_id="prot_1", protein_length=200, hits=[hit1, hit2])
        # Merged interval: 1..200 = full coverage
        assert pd.total_coverage == pytest.approx(1.0, abs=0.01)


# ── TSV output tests ────────────────────────────────────────────────────


class TestFormatTsv:
    def test_header_present(self) -> None:
        results = [ProteinDomains(protein_id="p1", protein_length=100, hits=[])]
        tsv = format_tsv(results)
        header = tsv.split("\n")[0]
        assert "protein_id" in header
        assert "domain_name" in header
        assert "evalue" in header

    def test_protein_with_no_domains(self) -> None:
        results = [ProteinDomains(protein_id="lonely", protein_length=500, hits=[])]
        tsv = format_tsv(results)
        lines = tsv.strip().split("\n")
        assert len(lines) == 2  # header + 1 data row
        assert "no domains found" in lines[1]

    def test_protein_with_domains(self) -> None:
        hit = DomainHit(
            protein_id="prot_1", protein_length=300,
            domain_name="PF00001", accession="PF00001.21",
            description="7tm_1", env_from=10, env_to=200,
            evalue=1.5e-30, score=100.5, db_name="Pfam-A",
        )
        results = [ProteinDomains(protein_id="prot_1", protein_length=300, hits=[hit])]
        tsv = format_tsv(results)
        lines = tsv.strip().split("\n")
        assert len(lines) == 2
        fields = lines[1].split("\t")
        assert fields[0] == "prot_1"
        assert fields[2] == "PF00001"
        assert "1.50e-30" in fields[8]

    def test_multiple_proteins_multiple_domains(self) -> None:
        hit1 = DomainHit(
            protein_id="p1", protein_length=300,
            domain_name="PF00001", accession="PF00001",
            description="dom1", env_from=1, env_to=100,
            evalue=1e-30, score=100.0, db_name="Pfam-A",
        )
        hit2 = DomainHit(
            protein_id="p1", protein_length=300,
            domain_name="PF00002", accession="PF00002",
            description="dom2", env_from=150, env_to=280,
            evalue=1e-20, score=80.0, db_name="Pfam-A",
        )
        results = [
            ProteinDomains(protein_id="p1", protein_length=300, hits=[hit1, hit2]),
            ProteinDomains(protein_id="p2", protein_length=100, hits=[]),
        ]
        tsv = format_tsv(results)
        lines = tsv.strip().split("\n")
        assert len(lines) == 4  # header + 2 hits for p1 + 1 "no domains" for p2


# ── Domain map visual tests ────────────────────────────────────────────


class TestFormatDomainMap:
    def test_no_domains(self) -> None:
        pd = ProteinDomains(protein_id="empty_prot", protein_length=300, hits=[])
        result = format_domain_map(pd, width=50)
        assert "empty_prot" in result
        assert "no domains found" in result
        assert "|" in result  # has boundary markers
        assert "-" in result  # has uncharacterized region

    def test_single_domain(self) -> None:
        hit = DomainHit(
            protein_id="prot_1", protein_length=200,
            domain_name="AAA", accession="PF00001",
            description="", env_from=50, env_to=150,
            evalue=1e-20, score=80.0, db_name="Pfam-A",
        )
        pd = ProteinDomains(protein_id="prot_1", protein_length=200, hits=[hit])
        result = format_domain_map(pd, width=60)
        assert "prot_1" in result
        assert "1 domain" in result
        assert "=" in result  # domain region
        assert "200" in result  # protein length at end

    def test_two_domains(self) -> None:
        hits = [
            DomainHit(
                protein_id="multi", protein_length=500,
                domain_name="SH2", accession="PF00017",
                description="", env_from=10, env_to=100,
                evalue=1e-30, score=100.0, db_name="Pfam-A",
            ),
            DomainHit(
                protein_id="multi", protein_length=500,
                domain_name="Kinase", accession="PF00069",
                description="", env_from=250, env_to=480,
                evalue=1e-40, score=120.0, db_name="Pfam-A",
            ),
        ]
        pd = ProteinDomains(protein_id="multi", protein_length=500, hits=hits)
        result = format_domain_map(pd, width=70)
        assert "2 domains" in result
        assert "=" in result
        assert "-" in result

    def test_zero_length_protein(self) -> None:
        pd = ProteinDomains(protein_id="zero", protein_length=0, hits=[])
        result = format_domain_map(pd, width=50)
        assert "zero" in result
        assert "0 aa" in result


class TestFormatAllMaps:
    def test_includes_header_and_legend(self) -> None:
        results = [ProteinDomains(protein_id="p1", protein_length=100, hits=[])]
        output = format_all_maps(results, width=50)
        assert "DOMAIN ARCHITECTURE MAP" in output
        assert "Legend" in output

    def test_multiple_proteins(self) -> None:
        results = [
            ProteinDomains(protein_id="alpha", protein_length=200, hits=[]),
            ProteinDomains(protein_id="beta", protein_length=300, hits=[]),
        ]
        output = format_all_maps(results, width=50)
        assert "alpha" in output
        assert "beta" in output
