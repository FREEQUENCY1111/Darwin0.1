"""
Stream — The reactive event backbone of the ecosphere.

Like water in a jar ecosystem, the stream carries nutrients
between organisms. No organism knows about any other —
they only know what nutrients they release and consume.

This is a pub/sub system disguised as nature.
"""

from __future__ import annotations

import asyncio
import logging
import time
from collections.abc import Awaitable, Callable
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

logger = logging.getLogger("darwin.water")


class NutrientType(str, Enum):
    """Types of nutrients that flow through the water."""

    # Sunlight hits the jar
    GENOME_LOADED = "genome.loaded"

    # Flora produce these
    GENES_CALLED = "genes.called"
    MINIGENES_FOUND = "minigenes.found"
    PROTEINS_FOUND = "proteins.found"
    TRNA_DETECTED = "trna.detected"
    RRNA_DETECTED = "rrna.detected"
    CRISPR_DETECTED = "crispr.detected"
    SIGNAL_PEPTIDES_FOUND = "signal_peptides.found"
    OPERONS_GROUPED = "operons.grouped"
    TAXONOMY_INFERRED = "taxonomy.inferred"
    PLASMIDS_CLASSIFIED = "plasmids.classified"
    MOBILE_ELEMENTS_FOUND = "mobile_elements.found"
    RESISTANCE_GENES_FOUND = "resistance_genes.found"
    PROPHAGES_DETECTED = "prophages.detected"
    BGC_DETECTED = "bgc.detected"

    # Microbiome produces these
    QC_COMPLETED = "qc.completed"
    ANNOTATION_READY = "annotation.ready"
    OUTPUT_WRITTEN = "output.written"

    # Distress signals (toxins in the water)
    ERROR = "error"
    WARNING = "warning"

    # Life signals
    HEARTBEAT = "heartbeat"
    EQUILIBRIUM = "equilibrium"


@dataclass
class Nutrient:
    """
    A single nutrient flowing through the water.

    Like dissolved minerals in a jar ecosystem — any organism
    that needs this nutrient will absorb it automatically.
    """

    type: NutrientType
    data: Any = None
    source: str = "unknown"
    timestamp: float = field(default_factory=time.time)
    correlation_id: str | None = None  # tracks a single run through the jar

    def __repr__(self) -> str:
        return f"Nutrient({self.type.value} from {self.source})"


# Type alias for nutrient consumers
Consumer = Callable[[Nutrient], Awaitable[None]]


class Stream:
    """
    The water stream — reactive event bus.

    Organisms register what they feed on. When those nutrients
    appear in the water, they automatically wake up and consume them.

    The stream supports:
    - Multiple consumers per nutrient type (many organisms eat the same thing)
    - Wildcard consumers (decomposers that feed on everything)
    - Nutrient history (the sediment layer — what flowed before)
    - Async by nature (organisms don't block each other)
    """

    def __init__(self) -> None:
        self._consumers: dict[NutrientType, list[Consumer]] = {}
        self._wildcard_consumers: list[Consumer] = []
        self._sediment: list[Nutrient] = []  # history of all nutrients
        self._lock = asyncio.Lock()
        self._equilibrium_event = asyncio.Event()
        self._nutrient_counts: dict[NutrientType, int] = {}

    def feeds_on(self, nutrient_type: NutrientType) -> Callable:
        """
        Decorator — declare what an organism feeds on.

        Usage:
            @stream.feeds_on(NutrientType.GENES_CALLED)
            async def consume_genes(nutrient: Nutrient):
                ...
        """

        def decorator(func: Consumer) -> Consumer:
            if nutrient_type not in self._consumers:
                self._consumers[nutrient_type] = []
            self._consumers[nutrient_type].append(func)
            logger.debug(f"🌿 {func.__qualname__} now feeds on {nutrient_type.value}")
            return func

        return decorator

    def feeds_on_everything(self) -> Callable:
        """Decorator for decomposers that monitor all nutrients."""

        def decorator(func: Consumer) -> Consumer:
            self._wildcard_consumers.append(func)
            logger.debug(f"🦠 {func.__qualname__} feeds on ALL nutrients")
            return func

        return decorator

    def subscribe(self, nutrient_type: NutrientType, consumer: Consumer) -> None:
        """Imperative subscription — organism starts feeding on a nutrient."""
        if nutrient_type not in self._consumers:
            self._consumers[nutrient_type] = []
        self._consumers[nutrient_type].append(consumer)

    def subscribe_all(self, consumer: Consumer) -> None:
        """Imperative wildcard subscription."""
        self._wildcard_consumers.append(consumer)

    async def release(self, nutrient: Nutrient) -> None:
        """
        Release a nutrient into the water.

        Every organism that feeds on this nutrient type will
        absorb it concurrently — like nutrients dissolving
        into water and being absorbed by all nearby roots.
        """
        async with self._lock:
            self._sediment.append(nutrient)
            self._nutrient_counts[nutrient.type] = self._nutrient_counts.get(nutrient.type, 0) + 1

        logger.info(f"💧 {nutrient.type.value} released by {nutrient.source}")

        # Check if equilibrium has been reached — this must happen
        # regardless of whether any consumers exist, because equilibrium
        # is a property of the water itself, not of the organisms.
        if nutrient.type == NutrientType.OUTPUT_WRITTEN:
            self._equilibrium_event.set()
            eq_nutrient = Nutrient(
                type=NutrientType.EQUILIBRIUM,
                data=self.get_sediment_summary(),
                source="water",
                correlation_id=nutrient.correlation_id,
            )
            self._sediment.append(eq_nutrient)
            self._nutrient_counts[NutrientType.EQUILIBRIUM] = 1

        # Gather all consumers for this nutrient
        consumers = list(self._consumers.get(nutrient.type, []))
        consumers.extend(self._wildcard_consumers)

        if not consumers:
            logger.debug(f"  ↳ No organisms feed on {nutrient.type.value}")
            return

        # All organisms absorb concurrently
        results = await asyncio.gather(
            *[consumer(nutrient) for consumer in consumers],
            return_exceptions=True,
        )

        # Check for toxins (errors)
        for i, result in enumerate(results):
            if isinstance(result, Exception):
                logger.error(
                    f"☠️ Toxin! {consumers[i].__qualname__} choked on "
                    f"{nutrient.type.value}: {result}"
                )
                # Release error as toxin into the water
                toxin = Nutrient(
                    type=NutrientType.ERROR,
                    data={"error": str(result), "nutrient": nutrient.type.value},
                    source=consumers[i].__qualname__,
                    correlation_id=nutrient.correlation_id,
                )
                self._sediment.append(toxin)

    async def wait_for_equilibrium(self, timeout: float = 300.0) -> bool:
        """Wait for the ecosystem to reach equilibrium (all work done)."""
        try:
            await asyncio.wait_for(self._equilibrium_event.wait(), timeout=timeout)
            return True
        except asyncio.TimeoutError:
            logger.warning("⏰ Ecosystem did not reach equilibrium in time")
            return False

    def get_sediment(self, nutrient_type: NutrientType | None = None) -> list[Nutrient]:
        """Read the sediment layer — history of all nutrients that flowed."""
        if nutrient_type:
            return [n for n in self._sediment if n.type == nutrient_type]
        return list(self._sediment)

    def get_sediment_summary(self) -> dict[str, int]:
        """How many of each nutrient flowed through."""
        return {k.value: v for k, v in self._nutrient_counts.items()}

    def reset(self) -> None:
        """Drain the water — fresh start."""
        self._sediment.clear()
        self._nutrient_counts.clear()
        self._equilibrium_event.clear()

    @property
    def is_at_equilibrium(self) -> bool:
        return self._equilibrium_event.is_set()
