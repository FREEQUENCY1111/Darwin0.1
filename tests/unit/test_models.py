"""Tests for core data models."""

from darwin.models import Contig, Feature, FeatureType, Genome, Strand


class TestFeature:
    def test_length(self):
        f = Feature(
            seq_id="c1", feature_type=FeatureType.CDS,
            start=100, end=399, strand=Strand.FORWARD,
        )
        assert f.length == 300

    def test_product_default(self):
        f = Feature(
            seq_id="c1", feature_type=FeatureType.CDS,
            start=1, end=100, strand=Strand.FORWARD,
        )
        assert f.product == "hypothetical protein"

    def test_product_set(self):
        f = Feature(
            seq_id="c1", feature_type=FeatureType.CDS,
            start=1, end=100, strand=Strand.FORWARD,
            attributes={"product": "DNA gyrase subunit A"},
        )
        assert f.product == "DNA gyrase subunit A"

    def test_locus_tag(self):
        f = Feature(
            seq_id="c1", feature_type=FeatureType.CDS,
            start=1, end=100, strand=Strand.FORWARD,
            attributes={"locus_tag": "DARWIN_00001"},
        )
        assert f.locus_tag == "DARWIN_00001"


class TestContig:
    def test_length(self):
        c = Contig(id="c1", sequence="ATGCATGC")
        assert c.length == 8

    def test_gc_content(self):
        c = Contig(id="c1", sequence="ATGCATGC")
        assert c.gc_content == 0.5

    def test_gc_content_empty(self):
        c = Contig(id="c1", sequence="")
        assert c.gc_content == 0.0


class TestGenome:
    def test_total_length(self, sample_genome: Genome):
        assert sample_genome.total_length > 0

    def test_num_contigs(self, sample_genome: Genome):
        assert sample_genome.num_contigs == 2

    def test_gc_content_range(self, sample_genome: Genome):
        gc = sample_genome.gc_content
        assert 0.0 < gc < 1.0

    def test_summary(self, sample_genome: Genome):
        s = sample_genome.summary()
        assert s["name"] == "test_genome"
        assert s["num_contigs"] == 2
        assert s["total_bp"] > 0
        assert "cds_count" in s

    def test_all_features_empty(self, sample_genome: Genome):
        assert sample_genome.all_features == []

    def test_all_features_populated(self, annotated_genome: Genome):
        assert len(annotated_genome.all_features) == 2
