"""
Flora — The producers of the ecosphere.

Like plants in a jar ecosystem, flora take in raw energy
(genome data) and produce complex organic molecules
(gene calls, protein annotations, RNA predictions).

Each plant feeds on specific nutrients from the water
and releases new nutrients when it grows.

Flora don't know about each other. They just grow.
"""

from darwin.flora.aragorn import AragornPlant
from darwin.flora.barrnap import BarrnapPlant
from darwin.flora.base import Organism
from darwin.flora.crispard import CRISPARd
from darwin.flora.isescan_plant import ISEScanPlant
from darwin.flora.minigene import MiniGeneHunter
from darwin.flora.mob_suite_plant import MobSuitePlant
from darwin.flora.operons import OperonGrouper
from darwin.flora.phylo_16s import PhyloIdentifier
from darwin.flora.prodigal import ProdigalPlant
from darwin.flora.pyhmmer_plant import PyhmmerPlant
from darwin.flora.signal_scanner import SignalScanner

__all__ = [
    "Organism",
    "ProdigalPlant",
    "PyhmmerPlant",
    "AragornPlant",
    "BarrnapPlant",
    "MiniGeneHunter",
    "CRISPARd",
    "SignalScanner",
    "OperonGrouper",
    "PhyloIdentifier",
    "MobSuitePlant",
    "ISEScanPlant",
]
