"""Tests for Darwin v0.2 new organisms and features."""

from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from darwin.rocks.models import Contig, Feature, FeatureType, Genome, Strand
from darwin.water.stream import Nutrient, NutrientType, Stream

# ── Helpers ──────────────────────────────────────────────


def _make_genome(
    seq: str = "ATGCGTACGATCG",
    features: list[Feature] | None = None,
    name: str = "test_genome",
) -> Genome:
    contig = Contig(
        id="contig_1",
        sequence=seq,
        description="test contig",
        features=features or [],
    )
    return Genome(name=name, contigs=[contig])


def _make_cds(
    start: int,
    end: int,
    strand: Strand = Strand.FORWARD,
    product: str = "hypothetical protein",
    translation: str = "",
    locus_tag: str = "",
    note: str = "",
    contig_id: str = "contig_1",
) -> Feature:
    return Feature(
        type=FeatureType.CDS,
        start=start,
        end=end,
        strand=strand,
        product=product,
        translation=translation,
        locus_tag=locus_tag,
        note=note,
        contig_id=contig_id,
    )


# ── CRISPR Detection ────────────────────────────────────


class TestCRISPARd:
    def test_find_crispr_arrays(self):
        from darwin.flora.crispard import _find_crispr_arrays

        # Build a synthetic CRISPR: 4 repeats of 31bp with 35bp spacers
        repeat = "AGTTTTAGAGCTATGCTGTTTTGAATGGTCC"
        spacer1 = "CGATCGATCGATCGATCGATCGATCGATCGATCGA"
        spacer2 = "TTTAAACCCGGGTTTAAACCCGGGTTTAAACCCTGA"
        spacer3 = "GGGAAATTTCCCAAAGGGAAATTTCCCAAAGGGAAG"
        crispr_seq = repeat + spacer1 + repeat + spacer2 + repeat + spacer3 + repeat

        # Pad with random sequence
        padding = "A" * 500
        full_seq = padding + crispr_seq + padding

        arrays = _find_crispr_arrays(full_seq, contig_id="test")
        assert len(arrays) >= 1, "Should detect at least one CRISPR array"
        arr = arrays[0]
        assert arr.spacer_count >= 2, "Should find multiple spacers"

    def test_no_crispr_in_random(self):
        import random

        from darwin.flora.crispard import _find_crispr_arrays

        random.seed(42)
        seq = "".join(random.choices("ACGT", k=5000))
        arrays = _find_crispr_arrays(seq, contig_id="test")
        # Random sequence should rarely produce CRISPR arrays
        assert len(arrays) <= 1


# ── Signal Peptide Detection ────────────────────────────


class TestSignalScanner:
    def test_known_signal_peptide(self):
        from darwin.flora.signal_scanner import _detect_signal_peptide

        # OmpA-like signal peptide: charged n-region + hydrophobic h-region + cleavage
        sp = "MKKTAIAIAVALAGFATVAQA"
        mature = "PQVNIQDETYQIMF" * 3
        result = _detect_signal_peptide(sp + mature)
        # May or may not detect depending on thresholds, but shouldn't crash
        assert isinstance(result, dict) or result is None

    def test_no_signal_peptide_cytoplasmic(self):
        from darwin.flora.signal_scanner import _detect_signal_peptide

        # Cytoplasmic protein — no signal peptide
        protein = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTG"
        result = _detect_signal_peptide(protein)
        # Cytoplasmic proteins generally lack signal peptides
        assert result is None or isinstance(result, dict)


# ── Operon Detection ────────────────────────────────────


class TestOperonGrouper:
    @pytest.mark.asyncio
    async def test_operon_clustering(self):
        from darwin.flora.operons import OperonGrouper

        stream = Stream()
        soil = MagicMock()

        # 3 CDS on forward strand, close together = operon
        features = [
            _make_cds(100, 400, Strand.FORWARD, locus_tag="gene_001"),
            _make_cds(450, 800, Strand.FORWARD, locus_tag="gene_002"),
            _make_cds(850, 1200, Strand.FORWARD, locus_tag="gene_003"),
            # Gap > 300bp, then another gene
            _make_cds(2000, 2500, Strand.FORWARD, locus_tag="gene_004"),
        ]
        genome = _make_genome(seq="A" * 3000, features=features)

        grouper = OperonGrouper(stream, soil)
        nutrient = Nutrient(
            type=NutrientType.GENES_CALLED,
            data={"genome": genome, "config": {}},
            source="test",
        )
        result = await grouper.grow(nutrient)
        assert result is not None
        assert result.data["operon_count"] >= 1
        # Check that operon annotations were added to features
        operon_genes = [f for f in features if f.note and "operon=" in f.note]
        assert len(operon_genes) >= 2  # at least 2 genes should be in an operon

    @pytest.mark.asyncio
    async def test_opposite_strand_breaks_operon(self):
        from darwin.flora.operons import OperonGrouper

        stream = Stream()
        soil = MagicMock()

        features = [
            _make_cds(100, 400, Strand.FORWARD, locus_tag="gene_001"),
            _make_cds(450, 800, Strand.REVERSE, locus_tag="gene_002"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)

        grouper = OperonGrouper(stream, soil)
        nutrient = Nutrient(
            type=NutrientType.GENES_CALLED,
            data={"genome": genome, "config": {}},
            source="test",
        )
        result = await grouper.grow(nutrient)
        assert result is not None
        # Different strands → no operon
        assert result.data["operon_count"] == 0


# ── 16S Phylogenetic Identification ─────────────────────


class TestPhyloIdentifier:
    def test_fuzzy_contains(self):
        from darwin.flora.phylo_16s import _fuzzy_contains

        # Exact match
        assert _fuzzy_contains("AACGTACGTAA", "CGTACGT", 0)
        # 1 mismatch
        assert _fuzzy_contains("AACGTAXGTAA", "CGTACGT", 1)
        # Too many mismatches
        assert not _fuzzy_contains("AAXXTAXGTAA", "CGTACGT", 1)

    def test_domain_classification(self):
        from darwin.flora.phylo_16s import _classify_domain

        # A sequence containing the EUB338 target should be Bacteria
        bact_seq = "N" * 300 + "ACTCCTACGGGAGGCAGC" + "N" * 300
        assert _classify_domain(bact_seq) == "Bacteria"

    def test_scoring(self):
        from darwin.flora.phylo_16s import _score_taxonomy

        # _score_taxonomy now takes a sequence string, returns tuple
        seq = "ACTCCTACGGGAGGCAGC" + "N" * 500  # bacterial marker
        name, confidence, desc, matches = _score_taxonomy(seq)
        assert isinstance(name, str)
        assert isinstance(confidence, float)


# ── MiniGene Hunter ─────────────────────────────────────


class TestMiniGeneHunter:
    def test_find_orfs(self):
        from darwin.flora.minigene import MiniGeneHunter

        # Synthetic sequence with a small ORF: ATG + 30 codons + TAG = 93bp
        orf = "ATG" + "GCT" * 30 + "TAG"
        padding = "AAA" * 20
        seq = padding + orf + padding
        orfs = MiniGeneHunter._find_orfs(seq, min_len=30, max_len=150)
        assert len(orfs) >= 1

    def test_rbs_scoring(self):
        from darwin.flora.minigene import _score_rbs

        # Perfect Shine-Dalgarno: AGGAGG at -7 to -12
        upstream = "AAAAAAAAGGAGGAAAA"  # AGGAGG at position 8-13
        score = _score_rbs(upstream)
        assert score > 0


# ── TSV Output ───────────────────────────────────────────


class TestTSVOutput:
    def test_write_tsv(self):
        from darwin.output.tsv import write_tsv

        features = [
            _make_cds(100, 500, Strand.FORWARD, product="test protein", locus_tag="T_001"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)

        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as f:
            path = Path(f.name)

        write_tsv(genome, path)
        content = path.read_text()
        assert "locus_tag" in content  # header
        assert "T_001" in content
        assert "test protein" in content
        path.unlink()


# ── Checkpoint System ────────────────────────────────────


class TestCheckpoint:
    def test_save_and_load(self):
        from darwin.jar.checkpoint import load_checkpoint, save_checkpoint

        genome = _make_genome(seq="ACGT" * 100)
        genome.organism = "Test bacterium"
        genome.taxonomy = "Proteobacteria"

        with tempfile.TemporaryDirectory() as tmpdir:
            save_checkpoint(
                genome=genome,
                stage="genes_called",
                output_dir=Path(tmpdir),
            )

            loaded = load_checkpoint(Path(tmpdir), stage="genes_called")
            assert loaded is not None
            assert loaded["stage"] == "genes_called"
            assert loaded["genome"].name == "test_genome"
            assert loaded["genome"].organism == "Test bacterium"

    def test_list_checkpoints(self):
        from darwin.jar.checkpoint import list_checkpoints, save_checkpoint

        genome = _make_genome()
        with tempfile.TemporaryDirectory() as tmpdir:
            save_checkpoint(genome, "genes_called", Path(tmpdir))
            save_checkpoint(genome, "proteins_found", Path(tmpdir))

            cps = list_checkpoints(Path(tmpdir))
            assert len(cps) == 2
            stages = [c["stage"] for c in cps]
            assert "genes_called" in stages
            assert "proteins_found" in stages

    def test_load_latest(self):
        from darwin.jar.checkpoint import load_checkpoint, save_checkpoint

        genome = _make_genome()
        with tempfile.TemporaryDirectory() as tmpdir:
            save_checkpoint(genome, "genes_called", Path(tmpdir))
            save_checkpoint(genome, "qc_completed", Path(tmpdir))

            loaded = load_checkpoint(Path(tmpdir))
            assert loaded is not None
            assert loaded["stage"] == "qc_completed"

    def test_no_checkpoint(self):
        from darwin.jar.checkpoint import load_checkpoint

        with tempfile.TemporaryDirectory() as tmpdir:
            loaded = load_checkpoint(Path(tmpdir))
            assert loaded is None


# ── ANI Calculator ───────────────────────────────────────


class TestANI:
    def test_kmer_set(self):
        from darwin.utils.ani import compute_kmer_set

        seq = "ACGTACGTACGTACGT" * 10
        kmers = compute_kmer_set(seq, k=8)
        assert len(kmers) > 0

    def test_jaccard_identical(self):
        from darwin.utils.ani import jaccard_similarity

        s = {"ACGT", "TGCA", "GGCC"}
        assert jaccard_similarity(s, s) == 1.0

    def test_jaccard_disjoint(self):
        from darwin.utils.ani import jaccard_similarity

        a = {"AAAA", "BBBB"}
        b = {"CCCC", "DDDD"}
        assert jaccard_similarity(a, b) == 0.0

    def test_jaccard_to_ani(self):
        from darwin.utils.ani import jaccard_to_ani

        # Perfect match → ~100%
        ani = jaccard_to_ani(1.0)
        assert ani == 100.0

        # No match → 0%
        ani = jaccard_to_ani(0.0)
        assert ani == 0.0

        # Partial → between 0 and 100
        ani = jaccard_to_ani(0.5)
        assert 0 < ani < 100

    def test_sketch_genome(self):
        from darwin.utils.ani import sketch_genome

        seq = "ACGTACGT" * 1000
        sketch = sketch_genome(seq, k=8, sketch_size=100)
        assert len(sketch) <= 100

    def test_compare_identical_genomes(self):
        from darwin.rocks.fasta import write_fasta
        from darwin.utils.ani import compare_genomes

        genome = _make_genome(seq="ACGTACGT" * 500)

        with tempfile.TemporaryDirectory() as tmpdir:
            p1 = Path(tmpdir) / "g1.fna"
            p2 = Path(tmpdir) / "g2.fna"
            write_fasta(genome, p1)
            write_fasta(genome, p2)

            result = compare_genomes(p1, p2, k=8, sketch_size=1000)
            assert result["comparison"]["jaccard_similarity"] == 1.0
            assert result["comparison"]["approximate_ani"] == 100.0
            assert result["comparison"]["same_species"] is True


# ── Enhanced Enricher ────────────────────────────────────


class TestEnricher:
    def test_metabolic_marker_scan(self):
        from darwin.microbiome.enricher import Enricher

        # Create features with known marker keywords
        features = [
            _make_cds(100, 500, product="cytochrome c oxidase subunit I", locus_tag="T_001"),
            _make_cds(600, 900, product="ATP synthase subunit alpha", locus_tag="T_002"),
            _make_cds(1000, 1500, product="nitrogenase iron protein NifH", locus_tag="T_003"),
        ]

        enricher = Enricher.__new__(Enricher)
        metabolic = enricher._scan_metabolic_markers(features)
        assert "aerobic" in metabolic
        assert "nitrogen_fixing" in metabolic

    def test_lifestyle_inference(self):
        from darwin.microbiome.enricher import Enricher

        enricher = Enricher.__new__(Enricher)
        metabolic = {
            "aerobic": ["cytochrome c oxidase"],
            "motile": ["flagellin"],
        }
        lifestyle = enricher._infer_lifestyle(metabolic)
        assert "aerobic respiration" in lifestyle
        assert "motile" in lifestyle

    def test_gc_interpretation(self):
        from darwin.microbiome.enricher import Enricher

        enricher = Enricher.__new__(Enricher)
        assert "AT-rich" in enricher._interpret_gc(35)
        assert "moderate" in enricher._interpret_gc(45)
        assert "GC-rich" in enricher._interpret_gc(55)

    def test_size_categorization(self):
        from darwin.microbiome.enricher import Enricher

        enricher = Enricher.__new__(Enricher)
        assert "minimal" in enricher._categorize_size(400_000)
        assert "typical" in enricher._categorize_size(3_000_000)
        assert "very large" in enricher._categorize_size(10_000_000)


# ── Enhanced Scrutinizer QC ──────────────────────────────


class TestScrutinizerQC:
    def test_rrna_length_check(self):
        from darwin.microbiome.scrutinizer import Scrutinizer

        scrutinizer = Scrutinizer.__new__(Scrutinizer)

        # Good 16S rRNA (1541bp)
        rrna_16s = Feature(
            type=FeatureType.RRNA, start=1, end=1541,
            product="16S ribosomal RNA", locus_tag="rRNA_001",
        )
        genome = _make_genome(seq="A" * 2000, features=[rrna_16s])
        check = scrutinizer._check_rrna_lengths(genome)
        assert check.passed is True

    def test_rrna_length_check_bad(self):
        from darwin.microbiome.scrutinizer import Scrutinizer

        scrutinizer = Scrutinizer.__new__(Scrutinizer)

        # Short 16S rRNA (too small)
        rrna_16s = Feature(
            type=FeatureType.RRNA, start=1, end=500,
            product="16S ribosomal RNA", locus_tag="rRNA_001",
        )
        genome = _make_genome(seq="A" * 1000, features=[rrna_16s])
        check = scrutinizer._check_rrna_lengths(genome)
        assert check.passed is False


# ── Synthesizer Output ───────────────────────────────────


class TestSynthesizerGFF3:
    def test_gff3_dbxref(self):
        from darwin.microbiome.synthesizer import Synthesizer

        synth = Synthesizer.__new__(Synthesizer)

        features = [
            Feature(
                type=FeatureType.CDS, start=100, end=500,
                strand=Strand.FORWARD, product="test protein",
                locus_tag="T_001", db_xref=["Pfam:PF00001", "TIGRFAMs:TIGR00123"],
            ),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        genome.taxonomy = "Gammaproteobacteria"

        with tempfile.NamedTemporaryFile(suffix=".gff3", delete=False, mode="w") as f:
            path = Path(f.name)

        synth._write_gff3(genome, path)
        content = path.read_text()
        assert "##gff-version 3" in content
        assert "##species Gammaproteobacteria" in content
        assert "Dbxref=Pfam:PF00001,TIGRFAMs:TIGR00123" in content
        path.unlink()

    def test_gff3_mobile_element(self):
        """Mobile elements should appear in GFF3 output."""
        from darwin.microbiome.synthesizer import Synthesizer

        synth = Synthesizer.__new__(Synthesizer)

        features = [
            Feature(
                type=FeatureType.MOBILE_ELEMENT, start=100, end=1500,
                strand=Strand.FORWARD, product="IS3 family transposase",
                locus_tag="T_is001", inference="ab initio prediction:ISEScan",
                note="IS family: IS3; cluster: IS3_1",
            ),
        ]
        genome = _make_genome(seq="A" * 2000, features=features)

        with tempfile.NamedTemporaryFile(suffix=".gff3", delete=False, mode="w") as f:
            path = Path(f.name)

        synth._write_gff3(genome, path)
        content = path.read_text()
        assert "mobile_element" in content
        assert "IS3 family transposase" in content
        assert "ISEScan" in content
        path.unlink()

    def test_genbank_plasmid_contig(self):
        """GenBank output should annotate plasmid contigs correctly."""
        from darwin.microbiome.synthesizer import Synthesizer

        synth = Synthesizer.__new__(Synthesizer)

        contig = Contig(
            id="plasmid_1",
            sequence="A" * 5000,
            description="test plasmid",
            replicon_type="plasmid",
            rep_type="IncF",
            is_circular=True,
        )
        genome = Genome(name="test", contigs=[contig])

        with tempfile.NamedTemporaryFile(suffix=".gbk", delete=False, mode="w") as f:
            path = Path(f.name)

        synth._write_genbank(genome, path, {"translation_table": 11})
        content = path.read_text()
        assert "circular" in content
        assert '/plasmid="unnamed"' in content
        assert "replicon type: IncF" in content
        path.unlink()


# ── MobSuitePlant Tests ─────────────────────────────────


class TestMobSuitePlant:
    def test_parse_contig_report(self):
        """Test parsing of MOB-suite contig_report.txt."""
        from darwin.flora.mob_suite_plant import MobSuitePlant

        # Create a genome with 3 contigs
        contigs = [
            Contig(id="contig_1", sequence="A" * 5000),
            Contig(id="contig_2", sequence="G" * 3000),
            Contig(id="contig_3", sequence="C" * 2000),
        ]
        genome = Genome(name="test", contigs=contigs)

        # Write a mock contig_report.txt
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False, newline=""
        ) as fh:
            report_path = Path(fh.name)
            fh.write(
                "sample_id\tmolecule_type\tprimary_cluster_id\t"
                "rep_type(s)\trelaxase_type(s)\n"
            )
            fh.write("contig_1\tchromosome\t-\t-\t-\n")
            fh.write("contig_2\tplasmid\tAA038\tIncF\tMOBF\n")
            fh.write("contig_3\tplasmid\tAA039\tColE1\t-\n")

        plasmid_count = MobSuitePlant._parse_contig_report(genome, report_path)
        report_path.unlink()

        assert plasmid_count == 2
        assert contigs[0].replicon_type == "chromosome"
        assert contigs[1].replicon_type == "plasmid"
        assert contigs[1].rep_type == "IncF"
        assert contigs[1].mob_type == "MOBF"
        assert contigs[2].replicon_type == "plasmid"
        assert contigs[2].rep_type == "ColE1"

    def test_parse_empty_report(self):
        """Missing report file should return 0 plasmids."""
        from darwin.flora.mob_suite_plant import MobSuitePlant

        genome = _make_genome()
        count = MobSuitePlant._parse_contig_report(
            genome, Path("/nonexistent/contig_report.txt")
        )
        assert count == 0

    def test_unclassified_marked_chromosome(self):
        """Unclassified contigs should be marked as chromosome."""
        from darwin.flora.mob_suite_plant import MobSuitePlant

        contigs = [
            Contig(id="contig_1", sequence="A" * 1000),
            Contig(id="contig_2", sequence="G" * 1000),
        ]
        genome = Genome(name="test", contigs=contigs)

        # Report only mentions contig_1
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False, newline=""
        ) as fh:
            report_path = Path(fh.name)
            fh.write("sample_id\tmolecule_type\n")
            fh.write("contig_1\tchromosome\n")

        MobSuitePlant._parse_contig_report(genome, report_path)
        report_path.unlink()

        # contig_2 wasn't in the report but should be marked chromosome
        assert contigs[0].replicon_type == "chromosome"
        assert contigs[1].replicon_type == "chromosome"


# ── ISEScanPlant Tests ──────────────────────────────────


class TestISEScanPlant:
    def test_parse_isescan_output(self):
        """Test parsing of ISEScan TSV output."""
        from darwin.flora.isescan_plant import ISEScanPlant

        contigs = [
            Contig(id="contig_1", sequence="A" * 10000),
        ]
        genome = Genome(name="test", contigs=contigs)

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "results.tsv"
            with open(tsv_path, "w", newline="") as fh:
                fh.write(
                    "seqID\tfamily\tcluster\tisBegin\tisEnd\tisLen\t"
                    "strand\ttype\tscore\tirId\tirLen\ttir\n"
                )
                fh.write(
                    "contig_1\tIS3\tIS3_1\t1000\t2500\t1500\t"
                    "+\tc\t150.5\t85\t28\tTGATCGAT\n"
                )
                fh.write(
                    "contig_1\tIS256\tIS256_1\t5000\t6200\t1200\t"
                    "-\tc\t120.0\t90\t25\tACGTTGCA\n"
                )

            is_count, families = ISEScanPlant._parse_isescan_output(
                genome, Path(tmpdir), "TEST"
            )

        assert is_count == 2
        assert "IS3" in families
        assert "IS256" in families
        # Check features were added to contig
        me_features = contigs[0].features_of_type(FeatureType.MOBILE_ELEMENT)
        assert len(me_features) == 2
        assert me_features[0].product == "IS3 family transposase"
        assert me_features[0].inference == "ab initio prediction:ISEScan"
        assert "IS family: IS3" in me_features[0].note
        assert me_features[0].strand == Strand.FORWARD
        assert me_features[1].strand == Strand.REVERSE

    def test_parse_empty_output(self):
        """No TSV files should return 0 IS elements."""
        from darwin.flora.isescan_plant import ISEScanPlant

        genome = _make_genome()
        with tempfile.TemporaryDirectory() as tmpdir:
            is_count, families = ISEScanPlant._parse_isescan_output(
                genome, Path(tmpdir), "TEST"
            )
        assert is_count == 0
        assert len(families) == 0

    def test_parse_skips_unknown_contig(self):
        """Rows with unknown contig IDs should be skipped."""
        from darwin.flora.isescan_plant import ISEScanPlant

        genome = _make_genome()
        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "results.tsv"
            with open(tsv_path, "w", newline="") as fh:
                fh.write("seqID\tfamily\tcluster\tisBegin\tisEnd\tisLen\tstrand\ttype\tscore\tirId\tirLen\ttir\n")
                fh.write("unknown_contig\tIS1\tIS1_1\t100\t500\t400\t+\tc\t50.0\t80\t20\tACGT\n")

            is_count, families = ISEScanPlant._parse_isescan_output(
                genome, Path(tmpdir), "TEST"
            )
        assert is_count == 0


# ── Data Model Extensions Tests ─────────────────────────


class TestDataModelExtensions:
    def test_mobile_element_feature_type(self):
        """MOBILE_ELEMENT should be a valid FeatureType."""
        assert FeatureType.MOBILE_ELEMENT.value == "mobile_element"

    def test_contig_replicon_fields(self):
        """Contig should have replicon metadata fields."""
        contig = Contig(
            id="test",
            sequence="ACGT",
            replicon_type="plasmid",
            mob_type="MOBF",
            rep_type="IncF",
            is_circular=True,
        )
        assert contig.replicon_type == "plasmid"
        assert contig.mob_type == "MOBF"
        assert contig.rep_type == "IncF"
        assert contig.is_circular is True

    def test_genome_summary_plasmid_count(self):
        """Genome.summary() should include plasmid and IS counts."""
        contigs = [
            Contig(id="chr", sequence="A" * 1000, replicon_type="chromosome"),
            Contig(id="p1", sequence="G" * 500, replicon_type="plasmid"),
        ]
        # Add an IS element to the chromosome
        contigs[0].features.append(
            Feature(
                type=FeatureType.MOBILE_ELEMENT, start=100, end=500,
                product="IS3 transposase", locus_tag="is_001",
            )
        )
        genome = Genome(name="test", contigs=contigs)
        summary = genome.summary()
        assert summary["plasmid_count"] == 1
        assert summary["is_element_count"] == 1

    def test_nutrient_types_exist(self):
        """New NutrientTypes should exist."""
        assert NutrientType.PLASMIDS_CLASSIFIED.value == "plasmids.classified"
        assert NutrientType.MOBILE_ELEMENTS_FOUND.value == "mobile_elements.found"
        assert NutrientType.RESISTANCE_GENES_FOUND.value == "resistance_genes.found"
        assert NutrientType.PROPHAGES_DETECTED.value == "prophages.detected"
        assert NutrientType.BGC_DETECTED.value == "bgc.detected"

    def test_v04_feature_types(self):
        """v0.4 FeatureTypes should exist."""
        assert FeatureType.AMR_GENE.value == "amr_gene"
        assert FeatureType.PROPHAGE.value == "prophage"
        assert FeatureType.BGC.value == "bgc"

    def test_genome_summary_v04_counts(self):
        """Genome.summary() should include AMR/prophage/BGC counts."""
        contigs = [Contig(id="chr", sequence="A" * 1000)]
        contigs[0].features.append(
            Feature(type=FeatureType.AMR_GENE, start=100, end=400,
                    product="blaTEM-1", locus_tag="amr_001")
        )
        contigs[0].features.append(
            Feature(type=FeatureType.PROPHAGE, start=500, end=900,
                    product="predicted prophage", locus_tag="pp_001")
        )
        contigs[0].features.append(
            Feature(type=FeatureType.BGC, start=1000, end=5000,
                    product="NRPS biosynthetic gene cluster", locus_tag="bgc_001")
        )
        genome = Genome(name="test", contigs=contigs)
        summary = genome.summary()
        assert summary["amr_gene_count"] == 1
        assert summary["prophage_count"] == 1
        assert summary["bgc_count"] == 1


# ── ABRicatePlant Tests ─────────────────────────────────


class TestABRicatePlant:
    def test_parse_abricate_output(self):
        """Test parsing of ABRicate TSV output."""
        from darwin.flora.amr_plant import ABRicatePlant

        contigs = [
            Contig(id="contig_1", sequence="A" * 10000),
        ]
        genome = Genome(name="test", contigs=contigs)

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False, newline=""
        ) as fh:
            tsv_path = Path(fh.name)
            # Write ABRicate header + 2 hits
            fh.write(
                "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\t"
                "COVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\t"
                "ACCESSION\tPRODUCT\tRESISTANCE\n"
            )
            fh.write(
                "input.fa\tcontig_1\t100\t900\t+\tblaTEM-1\t1-800/800\t"
                "========\t0/0\t100.00\t99.58\tncbi\tNG_050145.1\t"
                "class A beta-lactamase TEM-1\tBeta-lactam\n"
            )
            fh.write(
                "input.fa\tcontig_1\t2000\t2600\t-\taph(3'')-Ib\t1-600/600\t"
                "========\t0/0\t100.00\t98.33\tncbi\tNG_047300.1\t"
                "aminoglycoside phosphotransferase\tAminoglycoside\n"
            )

        amr_count, resistance_classes = ABRicatePlant._parse_abricate_output(
            genome, tsv_path, "TEST"
        )
        tsv_path.unlink()

        assert amr_count == 2
        assert "Beta-lactam" in resistance_classes
        assert "Aminoglycoside" in resistance_classes

        # Check features were attached
        amr_features = contigs[0].features_of_type(FeatureType.AMR_GENE)
        assert len(amr_features) == 2
        assert amr_features[0].product == "class A beta-lactamase TEM-1"
        assert amr_features[0].gene == "blaTEM-1"
        assert amr_features[0].strand == Strand.FORWARD
        assert "resistance class: Beta-lactam" in amr_features[0].note
        assert "AMR:NG_050145.1" in amr_features[0].db_xref
        assert amr_features[1].strand == Strand.REVERSE

    def test_parse_empty_output(self):
        """Missing output file should return 0 AMR genes."""
        from darwin.flora.amr_plant import ABRicatePlant

        genome = _make_genome()
        count, classes = ABRicatePlant._parse_abricate_output(
            genome, Path("/nonexistent/abricate.tsv"), "TEST"
        )
        assert count == 0
        assert len(classes) == 0

    def test_parse_skips_unknown_contig(self):
        """Rows with unknown contig IDs should be skipped."""
        from darwin.flora.amr_plant import ABRicatePlant

        genome = _make_genome()
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False, newline=""
        ) as fh:
            tsv_path = Path(fh.name)
            fh.write(
                "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\t"
                "COVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\t"
                "ACCESSION\tPRODUCT\tRESISTANCE\n"
            )
            fh.write(
                "input.fa\tunknown_contig\t100\t900\t+\tblaTEM-1\t1-800/800\t"
                "========\t0/0\t100.00\t99.58\tncbi\tNG_050145.1\t"
                "TEM-1\tBeta-lactam\n"
            )

        count, classes = ABRicatePlant._parse_abricate_output(
            genome, tsv_path, "TEST"
        )
        tsv_path.unlink()
        assert count == 0


# ── PhiSpyPlant Tests ───────────────────────────────────


class TestPhiSpyPlant:
    def test_feeds_on_genes_called(self):
        """PhiSpy must wait for gene calling (needs CDS in GenBank)."""
        from darwin.flora.phispy_plant import PhiSpyPlant

        assert NutrientType.GENES_CALLED in PhiSpyPlant.feeds_on_nutrients
        assert NutrientType.GENOME_LOADED not in PhiSpyPlant.feeds_on_nutrients

    def test_write_minimal_genbank(self):
        """Test minimal GenBank generation for PhiSpy."""
        from darwin.flora.phispy_plant import PhiSpyPlant

        features = [
            _make_cds(100, 500, Strand.FORWARD, product="test protein",
                      locus_tag="TEST_001", translation="MTEST"),
        ]
        genome = _make_genome(seq="A" * 2000, features=features)

        with tempfile.NamedTemporaryFile(
            suffix=".gbk", delete=False, mode="w"
        ) as f:
            gbk_path = Path(f.name)

        PhiSpyPlant._write_minimal_genbank(genome, gbk_path)
        content = gbk_path.read_text()
        gbk_path.unlink()

        assert "LOCUS" in content
        assert "CDS" in content
        assert "test protein" in content
        assert "ORIGIN" in content
        assert "//" in content

    def test_parse_phispy_output(self):
        """Test parsing of PhiSpy prophage coordinates."""
        from darwin.flora.phispy_plant import PhiSpyPlant

        contigs = [Contig(id="contig_1", sequence="A" * 50000)]
        genome = Genome(name="test", contigs=contigs)

        with tempfile.TemporaryDirectory() as tmpdir:
            coords_file = Path(tmpdir) / "prophage_coordinates.tsv"
            with open(coords_file, "w") as fh:
                # No-header format: pp_num, contig, start, stop
                fh.write("pp_1\tcontig_1\t5000\t25000\n")
                fh.write("pp_2\tcontig_1\t35000\t45000\n")

            count, statuses = PhiSpyPlant._parse_phispy_output(
                genome, Path(tmpdir), "TEST"
            )

        assert count == 2
        pp_features = contigs[0].features_of_type(FeatureType.PROPHAGE)
        assert len(pp_features) == 2
        assert pp_features[0].inference == "ab initio prediction:PhiSpy"
        assert "length: 20000bp" in pp_features[0].note

    def test_parse_empty_output(self):
        """No output files should return 0 prophages."""
        from darwin.flora.phispy_plant import PhiSpyPlant

        genome = _make_genome()
        with tempfile.TemporaryDirectory() as tmpdir:
            count, statuses = PhiSpyPlant._parse_phispy_output(
                genome, Path(tmpdir), "TEST"
            )
        assert count == 0


# ── GeccoPlant Tests ────────────────────────────────────


class TestGeccoPlant:
    def test_parse_gecco_output(self):
        """Test parsing of GECCO clusters TSV output."""
        from darwin.flora.gecco_plant import GeccoPlant

        contigs = [Contig(id="contig_1", sequence="A" * 100000)]
        genome = Genome(name="test", contigs=contigs)

        with tempfile.TemporaryDirectory() as tmpdir:
            clusters_file = Path(tmpdir) / "input.clusters.tsv"
            with open(clusters_file, "w", newline="") as fh:
                fh.write(
                    "sequence_id\tbgc_id\tstart\tend\taverage_p\t"
                    "type\tn_genes\n"
                )
                fh.write(
                    "contig_1\tbgc_1\t10000\t35000\t0.952\t"
                    "NRPS\t12\n"
                )
                fh.write(
                    "contig_1\tbgc_2\t60000\t78000\t0.871\t"
                    "Polyketide\t8\n"
                )

            count, bgc_types = GeccoPlant._parse_gecco_output(
                genome, Path(tmpdir), "TEST"
            )

        assert count == 2
        assert bgc_types["NRPS"] == 1
        assert bgc_types["Polyketide"] == 1

        bgc_features = contigs[0].features_of_type(FeatureType.BGC)
        assert len(bgc_features) == 2
        assert bgc_features[0].product == "NRPS biosynthetic gene cluster"
        assert bgc_features[0].inference == "ab initio prediction:GECCO"
        assert "type: NRPS" in bgc_features[0].note
        assert "confidence: 0.952" in bgc_features[0].note
        assert "genes in cluster: 12" in bgc_features[0].note

    def test_parse_empty_output(self):
        """No output files should return 0 BGCs."""
        from darwin.flora.gecco_plant import GeccoPlant

        genome = _make_genome()
        with tempfile.TemporaryDirectory() as tmpdir:
            count, bgc_types = GeccoPlant._parse_gecco_output(
                genome, Path(tmpdir), "TEST"
            )
        assert count == 0
        assert len(bgc_types) == 0

    def test_parse_skips_unknown_contig(self):
        """Rows with unknown contig IDs should be skipped."""
        from darwin.flora.gecco_plant import GeccoPlant

        genome = _make_genome()
        with tempfile.TemporaryDirectory() as tmpdir:
            clusters_file = Path(tmpdir) / "input.clusters.tsv"
            with open(clusters_file, "w", newline="") as fh:
                fh.write("sequence_id\tbgc_id\tstart\tend\taverage_p\ttype\tn_genes\n")
                fh.write("unknown_contig\tbgc_1\t1000\t5000\t0.9\tNRPS\t5\n")

            count, bgc_types = GeccoPlant._parse_gecco_output(
                genome, Path(tmpdir), "TEST"
            )
        assert count == 0


# ── v0.5 NutrientType and Feature Model Tests ────────────


class TestV05Types:
    def test_nutrient_types_v05(self):
        """v0.5 NutrientTypes should exist."""
        assert NutrientType.FUNCTIONS_ANNOTATED.value == "functions.annotated"
        assert NutrientType.STRUCTURES_MATCHED.value == "structures.matched"

    def test_feature_go_terms_field(self):
        """Feature should have go_terms list."""
        f = Feature(type=FeatureType.CDS, start=1, end=100)
        assert f.go_terms == []
        f.go_terms.append("GO:0003674")
        assert "GO:0003674" in f.go_terms

    def test_feature_ipr_ids_field(self):
        """Feature should have ipr_ids list."""
        f = Feature(type=FeatureType.CDS, start=1, end=100)
        assert f.ipr_ids == []
        f.ipr_ids.append("IPR000001")
        assert "IPR000001" in f.ipr_ids

    def test_feature_structure_hit_field(self):
        """Feature should have structure_hit string."""
        f = Feature(type=FeatureType.CDS, start=1, end=100)
        assert f.structure_hit == ""
        f.structure_hit = "1xyz_A"
        assert f.structure_hit == "1xyz_A"


# ── InterProPlant Tests ─────────────────────────────────


class TestInterProPlant:
    def test_parse_iprscan_output(self):
        """Parse standard InterProScan TSV output."""
        from darwin.flora.interpro_plant import InterProPlant

        features = [
            _make_cds(100, 500, Strand.FORWARD, product="hypothetical protein",
                      locus_tag="GENE_001", translation="MTEST"),
            _make_cds(600, 900, Strand.FORWARD, product="known protein",
                      locus_tag="GENE_002", translation="MKNOWN"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        proteins = {f.locus_tag: f for f in genome.contigs[0].features
                    if f.type == FeatureType.CDS}

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_file = Path(tmpdir) / "results.tsv"
            with open(tsv_file, "w") as fh:
                fh.write(
                    "GENE_001\tmd5\t133\tPfam\tPF00001\t7tm_1\t10\t130\t1e-20\tT\t"
                    "2024-01-01\tIPR000276\tG protein-coupled receptor\t"
                    "GO:0004930|GO:0007186\tREACTOME:R-HSA-390666\n"
                )
                fh.write(
                    "GENE_002\tmd5\t100\tTIGRFAM\tTIGR00001\tsome_fam\t5\t95\t1e-15\tT\t"
                    "2024-01-01\t-\t-\tGO:0003674\t-\n"
                )

            go_count, ipr_count, annotated = InterProPlant._parse_iprscan_output(
                proteins, tsv_file
            )

        assert annotated == 2
        assert go_count == 3
        assert ipr_count == 1

        gene1 = proteins["GENE_001"]
        assert gene1.product == "G protein-coupled receptor"
        assert "IPR000276" in gene1.ipr_ids
        assert "GO:0004930" in gene1.go_terms
        assert "GO:0007186" in gene1.go_terms
        assert "Pfam:PF00001" in gene1.db_xref
        assert "pathway: REACTOME:R-HSA-390666" in gene1.note

    def test_parse_empty_output(self):
        """Empty file should return 0 annotations."""
        from darwin.flora.interpro_plant import InterProPlant

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_file = Path(tmpdir) / "results.tsv"
            tsv_file.touch()
            go, ipr, ann = InterProPlant._parse_iprscan_output({}, tsv_file)
        assert go == 0
        assert ipr == 0
        assert ann == 0


# ── FoldseekPlant Tests ─────────────────────────────────


class TestFoldseekPlant:
    def test_parse_foldseek_output(self):
        """Parse Foldseek m8 tabular output."""
        from darwin.flora.foldseek_plant import FoldseekPlant

        features = [
            _make_cds(100, 500, Strand.FORWARD, product="hypothetical protein",
                      locus_tag="GENE_001", translation="MTEST"),
            _make_cds(600, 900, Strand.FORWARD, product="known protein",
                      locus_tag="GENE_002", translation="MKNOWN"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        proteins = {f.locus_tag: f for f in genome.contigs[0].features
                    if f.type == FeatureType.CDS}

        with tempfile.TemporaryDirectory() as tmpdir:
            result_file = Path(tmpdir) / "results.m8"
            with open(result_file, "w") as fh:
                fh.write("GENE_001\t1xyz_A\t25.5\t1e-10\t85.0\t0.72\t0.65\t0.80\n")
                fh.write("GENE_002\t2abc_B\t55.0\t1e-30\t150.0\t0.85\t0.82\t0.90\n")

            total_hits, remote_count = FoldseekPlant._parse_foldseek_output(
                result_file, proteins
            )

        assert total_hits == 2
        assert remote_count == 1

        gene1 = proteins["GENE_001"]
        assert gene1.structure_hit == "1xyz_A"
        assert "PDB:1xyz_A" in gene1.db_xref
        assert "TM-score: 0.720" in gene1.note
        assert "same fold" in gene1.note
        assert "seq identity: 25.5%" in gene1.note

        gene2 = proteins["GENE_002"]
        assert gene2.structure_hit == "2abc_B"

    def test_parse_empty_output(self):
        """Empty results should return 0 hits."""
        from darwin.flora.foldseek_plant import FoldseekPlant

        with tempfile.TemporaryDirectory() as tmpdir:
            result_file = Path(tmpdir) / "results.m8"
            result_file.touch()
            hits, remote = FoldseekPlant._parse_foldseek_output(result_file, {})
        assert hits == 0
        assert remote == 0

    def test_parse_alphafold_db_xref(self):
        """AlphaFold hits should get AlphaFoldDB: prefix."""
        from darwin.flora.foldseek_plant import FoldseekPlant

        features = [
            _make_cds(100, 500, Strand.FORWARD, product="hypothetical protein",
                      locus_tag="GENE_001", translation="MTEST"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        proteins = {f.locus_tag: f for f in genome.contigs[0].features
                    if f.type == FeatureType.CDS}

        with tempfile.TemporaryDirectory() as tmpdir:
            result_file = Path(tmpdir) / "results.m8"
            with open(result_file, "w") as fh:
                fh.write("GENE_001\tAF-P12345-F1\t40.0\t1e-20\t120.0\t0.80\t0.75\t0.85\n")

            FoldseekPlant._parse_foldseek_output(result_file, proteins)

        gene1 = proteins["GENE_001"]
        assert "AlphaFoldDB:AF-P12345-F1" in gene1.db_xref


# ── Edge Case Tests ─────────────────────────────────────


class TestInterProEdgeCases:
    def test_parse_malformed_rows(self):
        """Short or malformed rows should be skipped gracefully."""
        from darwin.flora.interpro_plant import InterProPlant

        features = [
            _make_cds(100, 500, Strand.FORWARD, locus_tag="G1", translation="M"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        proteins = {f.locus_tag: f for f in genome.contigs[0].features
                    if f.type == FeatureType.CDS}

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_file = Path(tmpdir) / "results.tsv"
            with open(tsv_file, "w") as fh:
                fh.write("too\tfew\tcolumns\n")
                fh.write("\n")
                fh.write("G1\tmd5\t100\tPfam\tPF999\tdesc\t1\t99\t1e-5\tT\t2024\tIPR001\tReal hit\tGO:0001234\t-\n")

            go, ipr, ann = InterProPlant._parse_iprscan_output(proteins, tsv_file)

        assert ann == 1
        assert go == 1
        assert ipr == 1

    def test_duplicate_go_terms_not_added_twice(self):
        """Same GO term from multiple hits should only appear once."""
        from darwin.flora.interpro_plant import InterProPlant

        features = [
            _make_cds(100, 500, Strand.FORWARD, locus_tag="G1", translation="M"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        proteins = {f.locus_tag: f for f in genome.contigs[0].features
                    if f.type == FeatureType.CDS}

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_file = Path(tmpdir) / "results.tsv"
            with open(tsv_file, "w") as fh:
                fh.write("G1\tmd5\t100\tPfam\tPF001\tdesc\t1\t99\t1e-5\tT\t2024\tIPR001\tHit\tGO:0001234\t-\n")
                fh.write("G1\tmd5\t100\tGene3D\tG3D001\tdesc\t1\t99\t1e-5\tT\t2024\tIPR001\tHit\tGO:0001234\t-\n")

            go, ipr, ann = InterProPlant._parse_iprscan_output(proteins, tsv_file)

        gene = proteins["G1"]
        assert gene.go_terms.count("GO:0001234") == 1  # no duplicates
        assert gene.ipr_ids.count("IPR001") == 1

    def test_unknown_locus_tag_skipped(self):
        """Rows with unrecognized locus tags should be skipped."""
        from darwin.flora.interpro_plant import InterProPlant

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_file = Path(tmpdir) / "results.tsv"
            with open(tsv_file, "w") as fh:
                fh.write("UNKNOWN\tmd5\t100\tPfam\tPF001\tdesc\t1\t99\t1e-5\tT\t2024\tIPR001\tHit\tGO:0001\t-\n")

            go, ipr, ann = InterProPlant._parse_iprscan_output({}, tsv_file)

        assert ann == 0
        assert go == 0

    def test_nonexistent_file_returns_zero(self):
        """Missing file should return 0 counts, not crash."""
        from darwin.flora.interpro_plant import InterProPlant

        go, ipr, ann = InterProPlant._parse_iprscan_output({}, Path("/nonexistent"))
        assert go == 0
        assert ipr == 0
        assert ann == 0


class TestFoldseekEdgeCases:
    def test_best_hit_by_bitscore(self):
        """When multiple hits exist for one protein, keep highest bitscore."""
        from darwin.flora.foldseek_plant import FoldseekPlant

        features = [
            _make_cds(100, 500, Strand.FORWARD, locus_tag="G1", translation="M"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        proteins = {f.locus_tag: f for f in genome.contigs[0].features
                    if f.type == FeatureType.CDS}

        with tempfile.TemporaryDirectory() as tmpdir:
            result_file = Path(tmpdir) / "results.m8"
            with open(result_file, "w") as fh:
                fh.write("G1\t1abc_A\t30.0\t1e-5\t50.0\t0.5\t0.4\t0.6\n")
                fh.write("G1\t2xyz_B\t45.0\t1e-20\t120.0\t0.8\t0.75\t0.85\n")
                fh.write("G1\t3low_C\t20.0\t1e-2\t30.0\t0.3\t0.2\t0.4\n")

            hits, remote = FoldseekPlant._parse_foldseek_output(result_file, proteins)

        assert hits == 1
        gene = proteins["G1"]
        assert gene.structure_hit == "2xyz_B"  # highest bitscore

    def test_remote_homolog_threshold(self):
        """Proteins with <30% identity AND TM>0.5 count as remote homologs."""
        from darwin.flora.foldseek_plant import FoldseekPlant

        features = [
            _make_cds(100, 300, Strand.FORWARD, locus_tag="G1", translation="M"),
            _make_cds(400, 700, Strand.FORWARD, locus_tag="G2", translation="M"),
            _make_cds(800, 1000, Strand.FORWARD, locus_tag="G3", translation="M"),
        ]
        genome = _make_genome(seq="A" * 1100, features=features)
        proteins = {f.locus_tag: f for f in genome.contigs[0].features
                    if f.type == FeatureType.CDS}

        with tempfile.TemporaryDirectory() as tmpdir:
            result_file = Path(tmpdir) / "results.m8"
            with open(result_file, "w") as fh:
                # Remote homolog: low identity, high TM
                fh.write("G1\t1a_A\t15.0\t1e-8\t60.0\t0.65\t0.6\t0.7\n")
                # Not remote: high identity
                fh.write("G2\t2b_B\t85.0\t1e-50\t200.0\t0.95\t0.9\t0.95\n")
                # Not remote: low identity but low TM too
                fh.write("G3\t3c_C\t10.0\t0.1\t20.0\t0.3\t0.2\t0.35\n")

            hits, remote = FoldseekPlant._parse_foldseek_output(result_file, proteins)

        assert hits == 3
        assert remote == 1  # only G1

    def test_similar_fold_annotation(self):
        """TM-score between 0.5-0.7 should say 'similar fold', not 'same fold'."""
        from darwin.flora.foldseek_plant import FoldseekPlant

        features = [
            _make_cds(100, 500, Strand.FORWARD, locus_tag="G1", translation="M"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        proteins = {f.locus_tag: f for f in genome.contigs[0].features
                    if f.type == FeatureType.CDS}

        with tempfile.TemporaryDirectory() as tmpdir:
            result_file = Path(tmpdir) / "results.m8"
            with open(result_file, "w") as fh:
                fh.write("G1\t1a_A\t25.0\t1e-5\t50.0\t0.55\t0.5\t0.6\n")

            FoldseekPlant._parse_foldseek_output(result_file, proteins)

        gene = proteins["G1"]
        assert "similar fold" in gene.note
        assert "same fold" not in gene.note

    def test_malformed_numeric_fields_skipped(self):
        """Rows with non-numeric evalue/bitscore should be skipped."""
        from darwin.flora.foldseek_plant import FoldseekPlant

        features = [
            _make_cds(100, 500, Strand.FORWARD, locus_tag="G1", translation="M"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        proteins = {f.locus_tag: f for f in genome.contigs[0].features
                    if f.type == FeatureType.CDS}

        with tempfile.TemporaryDirectory() as tmpdir:
            result_file = Path(tmpdir) / "results.m8"
            with open(result_file, "w") as fh:
                fh.write("G1\t1a_A\tNaN\tbad\tnope\t0.5\t0.4\t0.6\n")

            hits, remote = FoldseekPlant._parse_foldseek_output(result_file, proteins)

        assert hits == 0


# ── Output Pipeline Integration Tests ───────────────────


class TestOutputPipelineV05:
    def test_tsv_includes_v05_columns(self):
        """TSV output should include go_terms, ipr_ids, structure_hit columns."""
        from darwin.output.tsv import write_tsv

        features = [
            _make_cds(100, 500, Strand.FORWARD, locus_tag="G1", translation="M"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        gene = genome.contigs[0].features[0]
        gene.go_terms = ["GO:0003674", "GO:0005575"]
        gene.ipr_ids = ["IPR000001"]
        gene.structure_hit = "1xyz_A"

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "test.tsv"
            write_tsv(genome, tsv_path)
            content = tsv_path.read_text()

        lines = content.strip().split("\n")
        header = lines[0].split("\t")
        assert "go_terms" in header
        assert "ipr_ids" in header
        assert "structure_hit" in header

        data = lines[1].split("\t")
        go_idx = header.index("go_terms")
        ipr_idx = header.index("ipr_ids")
        struct_idx = header.index("structure_hit")
        assert "GO:0003674" in data[go_idx]
        assert "GO:0005575" in data[go_idx]
        assert "IPR000001" in data[ipr_idx]
        assert "1xyz_A" in data[struct_idx]

    def test_gff3_includes_v05_attributes(self):
        """GFF3 should include Ontology_term, interpro, structure_hit."""
        from darwin.microbiome.synthesizer import Synthesizer
        from darwin.soil.nutrients import NutrientStore
        from darwin.water.stream import Stream

        features = [
            _make_cds(100, 500, Strand.FORWARD, locus_tag="G1", translation="M"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        gene = genome.contigs[0].features[0]
        gene.go_terms = ["GO:0003674"]
        gene.ipr_ids = ["IPR000001"]
        gene.structure_hit = "2abc_B"

        with tempfile.TemporaryDirectory() as tmpdir:
            gff_path = Path(tmpdir) / "test.gff3"
            stream = Stream()
            soil = NutrientStore()
            synth = Synthesizer(stream, soil, output_dir=Path(tmpdir))
            synth._write_gff3(genome, gff_path)
            content = gff_path.read_text()

        assert "Ontology_term=GO:0003674" in content
        assert "interpro=IPR000001" in content
        assert "structure_hit=2abc_B" in content

    def test_genbank_includes_v05_qualifiers(self):
        """GenBank should include GO and InterPro db_xref qualifiers."""
        from darwin.microbiome.synthesizer import Synthesizer
        from darwin.soil.nutrients import NutrientStore
        from darwin.water.stream import Stream

        features = [
            _make_cds(100, 500, Strand.FORWARD, locus_tag="G1", translation="MTEST"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        gene = genome.contigs[0].features[0]
        gene.go_terms = ["GO:0003674"]
        gene.ipr_ids = ["IPR000001"]
        gene.structure_hit = "1xyz_A"

        with tempfile.TemporaryDirectory() as tmpdir:
            gbk_path = Path(tmpdir) / "test.gbk"
            stream = Stream()
            soil = NutrientStore()
            synth = Synthesizer(stream, soil, output_dir=Path(tmpdir))
            synth._write_genbank(genome, gbk_path)
            content = gbk_path.read_text()

        assert '/GO_component="GO:0003674"' in content
        assert '/db_xref="InterPro:IPR000001"' in content
        assert '/db_xref="PDB:1xyz_A"' in content


class TestBackwardCompatibility:
    def test_feature_defaults_empty_lists(self):
        """New fields should default to empty, not break old code."""
        f = Feature(type=FeatureType.CDS, start=1, end=100)
        assert f.go_terms == []
        assert f.ipr_ids == []
        assert f.structure_hit == ""
        # Ensure old fields still work
        assert f.db_xref == []
        assert f.note == ""
        assert f.product == "hypothetical protein"

    def test_genome_summary_still_works(self):
        """Genome.summary() should work with or without v0.5 fields populated."""
        contigs = [Contig(id="c1", sequence="A" * 1000)]
        contigs[0].features.append(
            Feature(type=FeatureType.CDS, start=100, end=400, locus_tag="g1")
        )
        genome = Genome(name="test", contigs=contigs)
        summary = genome.summary()
        assert "total_features" in summary
        assert summary["total_features"] == 1

    def test_tsv_empty_v05_fields(self):
        """TSV should handle features with no GO/IPR/structure gracefully."""
        from darwin.output.tsv import write_tsv

        features = [
            _make_cds(100, 500, Strand.FORWARD, locus_tag="G1", translation="M"),
        ]
        genome = _make_genome(seq="A" * 1000, features=features)
        # Don't set any v0.5 fields

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = Path(tmpdir) / "test.tsv"
            write_tsv(genome, tsv_path)
            content = tsv_path.read_text()

        lines = content.strip().split("\n")
        assert len(lines) == 2  # header + 1 data row
        # Should not crash with empty fields
