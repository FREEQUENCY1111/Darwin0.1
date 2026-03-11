"""Annotator modules — each wraps an external bioinformatics tool."""

from darwin.annotators.prodigal import ProdigalAnnotator
from darwin.annotators.pyhmmer_annotator import PyhmmerAnnotator
from darwin.annotators.aragorn import AragornAnnotator
from darwin.annotators.barrnap import BarrnapAnnotator

__all__ = [
    "ProdigalAnnotator",
    "PyhmmerAnnotator",
    "AragornAnnotator",
    "BarrnapAnnotator",
]
