"""
Water — The reactive event stream.

Nothing in the ecosphere communicates directly.
Everything flows through the water. Nutrients, signals,
waste products — all carried by the current.

Components don't call each other. They:
  1. Release nutrients into the water (emit events)
  2. Feed on nutrients they need (subscribe to events)

The water is the ONLY connection between living things.
"""

from darwin.water.cycle import WaterCycle
from darwin.water.stream import Nutrient, NutrientType, Stream

__all__ = ["Stream", "Nutrient", "NutrientType", "WaterCycle"]
