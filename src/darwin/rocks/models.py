"""
Core data models — the bedrock of the ecosphere.

These models are intentionally simple and immutable once created.
Like rocks in a jar — they provide structure but don't act.
"""

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
    TMRNA = "tmRNA"
    RRNA = "rRNA"
    REPEAT = "repeat_region"
    CRISPR = "CRISPR"
    SIGNAL_PEPTIDE = "signal_peptide"
    MOBILE_ELEMENT = "mobile_element"
    MISC = "misc_feature"


@dataclass
class Feature:
    """A single genomic feature — a pebble on the rock."""

    type: FeatureType
    start: int
    end: int
    strand: Strand = Strand.UNKNOWN
    score: float = 0.0
    contig_id: str = ""
    locus_tag: str = ""
    product: str = "hypothetical protein"
    inference: str = ""
    translation: str = ""
    gene: str = ""
    note: str = ""
    db_xref: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return abs(self.end - self.start)

    @property
    def location_str(self) -> str:
        if self.strand == Strand.REVERSE:
            return f"complement({self.start}..{self.end})"
        return f"{self.start}..{self.end}"


@dataclass
class Contig:
    """A single contiguous sequence — a layer of rock."""

    id: str
    sequence: str
    description: str = ""
    features: list[Feature] = field(default_factory=list)

    # Replicon classification (populated by MobSuitePlant)
    replicon_type: str = ""        # "chromosome" | "plasmid" | "" (unclassified)
    mob_type: str = ""             # MOB-suite mobility class (MOBP, MOBF, etc.)
    rep_type: str = ""             # Replicon incompatibility type (IncF, ColE1, etc.)
    is_circular: bool | None = None  # None = unknown

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        if not self.sequence:
            return 0.0
        gc = sum(1 for b in self.sequence.upper() if b in "GC")
        return round(gc / len(self.sequence) * 100, 2)

    def features_of_type(self, ftype: FeatureType) -> list[Feature]:
        return [f for f in self.features if f.type == ftype]


@dataclass
class Genome:
    """
    The whole rock — the complete genome.

    This is what sunlight hits when it enters the jar.
    Everything in the ecosphere ultimately references this.
    """

    name: str
    contigs: list[Contig] = field(default_factory=list)
    source_file: Path | None = None
    organism: str = ""
    taxonomy: str = ""

    @property
    def total_length(self) -> int:
        return sum(c.length for c in self.contigs)

    @property
    def num_contigs(self) -> int:
        return len(self.contigs)

    @property
    def gc_content(self) -> float:
        total_gc = sum(sum(1 for b in c.sequence.upper() if b in "GC") for c in self.contigs)
        total_len = self.total_length
        if total_len == 0:
            return 0.0
        return round(total_gc / total_len * 100, 2)

    @property
    def all_features(self) -> list[Feature]:
        return [f for c in self.contigs for f in c.features]

    def features_of_type(self, ftype: FeatureType) -> list[Feature]:
        return [f for f in self.all_features if f.type == ftype]

    def summary(self) -> dict:
        features = self.all_features
        by_type: dict[str, int] = {}
        for f in features:
            by_type[f.type.value] = by_type.get(f.type.value, 0) + 1
        plasmid_count = sum(1 for c in self.contigs if c.replicon_type == "plasmid")
        is_count = by_type.get(FeatureType.MOBILE_ELEMENT.value, 0)
        return {
            "name": self.name,
            "organism": self.organism,
            "total_length_bp": self.total_length,
            "num_contigs": self.num_contigs,
            "gc_content": self.gc_content,
            "total_features": len(features),
            "features_by_type": by_type,
            "plasmid_count": plasmid_count,
            "is_element_count": is_count,
        }


@dataclass
class AnnotationConfig:
    """Configuration for an annotation run — like the recipe for filling the jar."""

    input_file: Path
    output_dir: Path = Path("darwin_output")
    locus_tag_prefix: str = "DARWIN"
    translation_table: int = 11  # bacterial
    evalue_threshold: float = 1e-10
    min_contig_length: int = 200
    cpus: int = 1
    hmm_databases: list[Path] = field(default_factory=list)
    skip_tools: list[str] = field(default_factory=list)
    metagenome_mode: bool = False

    def __post_init__(self) -> None:
        self.output_dir = Path(self.output_dir)
        self.input_file = Path(self.input_file)
