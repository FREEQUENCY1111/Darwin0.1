"""
Microbiome — The decomposers and refiners.

Like bacteria in a jar ecosystem that break down dead matter
into nutrients, the microbiome takes raw annotation output
and refines it into meaningful, quality-checked results.

They don't produce new data — they transform and validate.
"""

from darwin.microbiome.scrutinizer import Scrutinizer
from darwin.microbiome.enricher import Enricher
from darwin.microbiome.synthesizer import Synthesizer

__all__ = ["Scrutinizer", "Enricher", "Synthesizer"]
