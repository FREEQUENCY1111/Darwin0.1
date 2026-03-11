"""
Organism — The base class for all living things in the ecosphere.

Every organism has:
  - A name
  - What it feeds on (nutrient types it consumes)
  - What it produces (nutrient types it releases)
  - A reference to the water (stream) it lives in
  - A reference to the soil it draws from

Organisms don't control anything. They react.
When their food appears in the water, they eat and grow.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import Optional

from darwin.water.stream import Stream, Nutrient, NutrientType
from darwin.soil.nutrients import NutrientStore
from darwin.rocks.models import Genome

logger = logging.getLogger("darwin.flora")


class Organism(ABC):
    """
    Base class for all living things in the jar.

    Subclasses define:
    - feeds_on: what nutrient types wake this organism up
    - produces: what nutrient types this organism releases
    - grow(): the actual work this organism does
    """

    name: str = "unnamed"
    feeds_on_nutrients: list[NutrientType] = []
    produces_nutrients: list[NutrientType] = []

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        self.stream = stream
        self.soil = soil
        self.logger = logging.getLogger(f"darwin.flora.{self.name}")
        self._alive = True

    def plant(self) -> None:
        """
        Plant this organism in the jar.

        Registers with the water stream to receive nutrients
        it feeds on. Like roots extending into the water.
        """
        for nutrient_type in self.feeds_on_nutrients:
            self.stream.subscribe(nutrient_type, self._consume)
            self.logger.debug(
                f"🌱 {self.name} planted, feeding on {nutrient_type.value}"
            )

    async def _consume(self, nutrient: Nutrient) -> None:
        """
        Called when food appears in the water.

        The organism absorbs the nutrient and grows.
        """
        if not self._alive:
            return

        self.logger.info(f"🌿 {self.name} absorbing {nutrient.type.value}")

        try:
            result = await self.grow(nutrient)
            if result:
                await self.stream.release(result)
        except Exception as e:
            self.logger.error(f"🥀 {self.name} wilted: {e}")
            await self.stream.release(Nutrient(
                type=NutrientType.ERROR,
                data={"organism": self.name, "error": str(e)},
                source=self.name,
                correlation_id=nutrient.correlation_id,
            ))

    @abstractmethod
    async def grow(self, nutrient: Nutrient) -> Optional[Nutrient]:
        """
        Do the work. Absorb input, produce output.

        Returns a new Nutrient to release into the water,
        or None if nothing to produce.
        """
        ...

    def can_grow(self) -> bool:
        """Check if this organism has what it needs in the soil."""
        return True

    def wilt(self) -> None:
        """This organism dies — stops consuming."""
        self._alive = False
        self.logger.info(f"🥀 {self.name} has wilted")
