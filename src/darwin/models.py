"""Core data models shared across all Darwin modules."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


class Strand(str, Enum):
    FORWARD = "+"
    REVERSE = "-"
    UNKNOWN = "."


class FeatureType(str, Enum):
    CDS = "CDS"
    TRNA = "tRNA"
    RRNA = "rRNA"
    TMRNA = "tmRNA"
    REPEAT = "repeat_region"
    MISC = "misc_feature"


@dataclass
class Feature:
    """A single genomic feature (gene, tRNA, rRNA, etc.)."""

    seq_id: str
    feature_type: FeatureType
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    strand: Strand
    score: float | None = None
    phase: int | None = None  # reading frame for CDS (0, 1, 2)
    attributes: dict[str, str] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return abs(self.end - self.start) + 1

    @property
    def locus_tag(self) -> str:
        return self.attributes.get("locus_tag", "")

    @property
    def product(self) -> str:
        return self.attributes.get("product", "hypothetical protein")


@dataclass
class Contig:
    """A single contig/sequence from the input genome."""

    id: str
    sequence: str
    description: str = ""
    features: list[Feature] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        if not self.sequence:
            return 0.0
        gc = sum(1 for b in self.sequence.upper() if b in ("G", "C"))
        return gc / len(self.sequence)


@dataclass
class Genome:
    """Complete genome assembly with all contigs and annotations."""

    name: str
    contigs: list[Contig] = field(default_factory=list)
    metadata: dict[str, str] = field(default_factory=dict)

    @property
    def total_length(self) -> int:
        return sum(c.length for c in self.contigs)

    @property
    def num_contigs(self) -> int:
        return len(self.contigs)

    @property
    def all_features(self) -> list[Feature]:
        return [f for c in self.contigs for f in c.features]

    @property
    def gc_content(self) -> float:
        total_seq = "".join(c.sequence for c in self.contigs)
        if not total_seq:
            return 0.0
        gc = sum(1 for b in total_seq.upper() if b in ("G", "C"))
        return gc / len(total_seq)

    def summary(self) -> dict:
        features = self.all_features
        return {
            "name": self.name,
            "total_bp": self.total_length,
            "num_contigs": self.num_contigs,
            "gc_content": round(self.gc_content, 4),
            "total_features": len(features),
            "cds_count": sum(1 for f in features if f.feature_type == FeatureType.CDS),
            "trna_count": sum(1 for f in features if f.feature_type == FeatureType.TRNA),
            "rrna_count": sum(1 for f in features if f.feature_type == FeatureType.RRNA),
            "tmrna_count": sum(1 for f in features if f.feature_type == FeatureType.TMRNA),
        }


@dataclass
class AnnotationConfig:
    """Configuration for a Darwin annotation run."""

    input_path: Path
    output_dir: Path
    locus_tag_prefix: str = "DARWIN"
    translation_table: int = 11  # bacterial genetic code
    evalue: float = 1e-6
    cpus: int = 1
    metagenome: bool = False
    gram: str = ""  # "", "pos", "neg" for SignalP
    kingdom: str = "bac"  # "bac" or "arc"
    min_contig_len: int = 200
    formats: list[str] = field(default_factory=lambda: ["gff3", "gbk", "faa", "fna", "json"])
