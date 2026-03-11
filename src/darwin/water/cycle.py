"""
WaterCycle — Tracks the flow of nutrients through the ecosystem.

Like tracing dye through water to see where nutrients go.
This provides observability without control.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

from darwin.water.stream import Nutrient, NutrientType


@dataclass
class FlowRecord:
    """A single step in the nutrient flow."""

    nutrient_type: str
    source: str
    timestamp: float
    correlation_id: str | None = None


@dataclass
class WaterCycle:
    """
    Tracks the full lifecycle of a single run through the ecosphere.

    Not a controller — just an observer. Like measuring
    dissolved oxygen levels without disturbing the water.
    """

    correlation_id: str
    started_at: float = field(default_factory=time.time)
    ended_at: float | None = None
    flow: list[FlowRecord] = field(default_factory=list)
    errors: list[dict] = field(default_factory=list)
    reached_equilibrium: bool = False

    def record(self, nutrient: Nutrient) -> None:
        """Record a nutrient flowing past."""
        self.flow.append(
            FlowRecord(
                nutrient_type=nutrient.type.value,
                source=nutrient.source,
                timestamp=nutrient.timestamp,
                correlation_id=nutrient.correlation_id,
            )
        )
        if nutrient.type == NutrientType.ERROR:
            self.errors.append(nutrient.data)
        if nutrient.type == NutrientType.EQUILIBRIUM:
            self.reached_equilibrium = True
            self.ended_at = time.time()

    @property
    def duration(self) -> float | None:
        if self.ended_at:
            return round(self.ended_at - self.started_at, 2)
        return None

    @property
    def nutrients_flowed(self) -> list[str]:
        return [r.nutrient_type for r in self.flow]

    def summary(self) -> dict:
        return {
            "correlation_id": self.correlation_id,
            "duration_seconds": self.duration,
            "equilibrium": self.reached_equilibrium,
            "nutrients_flowed": len(self.flow),
            "errors": len(self.errors),
            "flow": [{"type": r.nutrient_type, "source": r.source} for r in self.flow],
        }
