"""Tests for the Water layer — reactive event stream."""

import pytest

from darwin.water.cycle import WaterCycle
from darwin.water.stream import Nutrient, NutrientType


class TestStream:
    @pytest.mark.asyncio
    async def test_release_and_consume(self, stream):
        received = []

        async def consumer(nutrient: Nutrient):
            received.append(nutrient)

        stream.subscribe(NutrientType.GENOME_LOADED, consumer)

        await stream.release(
            Nutrient(
                type=NutrientType.GENOME_LOADED,
                data={"test": True},
                source="test",
            )
        )

        assert len(received) == 1
        assert received[0].data == {"test": True}

    @pytest.mark.asyncio
    async def test_multiple_consumers(self, stream):
        results = []

        async def consumer_a(n: Nutrient):
            results.append("a")

        async def consumer_b(n: Nutrient):
            results.append("b")

        stream.subscribe(NutrientType.GENES_CALLED, consumer_a)
        stream.subscribe(NutrientType.GENES_CALLED, consumer_b)

        await stream.release(
            Nutrient(
                type=NutrientType.GENES_CALLED,
                source="test",
            )
        )

        assert "a" in results
        assert "b" in results

    @pytest.mark.asyncio
    async def test_wildcard_consumer(self, stream):
        all_nutrients = []

        async def monitor(n: Nutrient):
            all_nutrients.append(n.type)

        stream.subscribe_all(monitor)

        await stream.release(Nutrient(type=NutrientType.GENOME_LOADED, source="t"))
        await stream.release(Nutrient(type=NutrientType.GENES_CALLED, source="t"))

        assert NutrientType.GENOME_LOADED in all_nutrients
        assert NutrientType.GENES_CALLED in all_nutrients

    @pytest.mark.asyncio
    async def test_sediment(self, stream):
        await stream.release(Nutrient(type=NutrientType.GENOME_LOADED, source="t"))
        await stream.release(Nutrient(type=NutrientType.GENES_CALLED, source="t"))

        sediment = stream.get_sediment()
        assert len(sediment) >= 2

        genome_sediment = stream.get_sediment(NutrientType.GENOME_LOADED)
        assert len(genome_sediment) == 1

    @pytest.mark.asyncio
    async def test_equilibrium_on_output_written(self, stream):
        assert not stream.is_at_equilibrium

        await stream.release(
            Nutrient(
                type=NutrientType.OUTPUT_WRITTEN,
                source="test",
            )
        )

        assert stream.is_at_equilibrium

    def test_reset(self, stream):
        stream.reset()
        assert not stream.is_at_equilibrium
        assert stream.get_sediment() == []


class TestWaterCycle:
    def test_record(self):
        cycle = WaterCycle(correlation_id="test-123")
        nutrient = Nutrient(type=NutrientType.GENOME_LOADED, source="test")
        cycle.record(nutrient)
        assert len(cycle.flow) == 1
        assert cycle.nutrients_flowed == ["genome.loaded"]

    def test_summary(self):
        cycle = WaterCycle(correlation_id="test-456")
        cycle.record(Nutrient(type=NutrientType.GENOME_LOADED, source="a"))
        cycle.record(Nutrient(type=NutrientType.GENES_CALLED, source="b"))

        s = cycle.summary()
        assert s["correlation_id"] == "test-456"
        assert s["nutrients_flowed"] == 2
