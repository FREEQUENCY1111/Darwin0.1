"""
Rocks — The substrate. The immutable foundation.

Rocks don't move. Rocks don't change. Rocks don't care.
Everything else in the jar rests on them.

These are the core data models — Genome, Contig, Feature.
They are what they are. Organisms build on top of them
but never alter them.
"""

from darwin.rocks.fasta import parse_fasta, write_fasta, write_proteins
from darwin.rocks.models import (
    AnnotationConfig,
    Contig,
    Feature,
    FeatureType,
    Genome,
    Strand,
)

__all__ = [
    "Genome",
    "Contig",
    "Feature",
    "Strand",
    "FeatureType",
    "AnnotationConfig",
    "parse_fasta",
    "write_fasta",
    "write_proteins",
]
