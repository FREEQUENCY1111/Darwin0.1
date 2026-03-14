"""
Ecosphere — The mason jar that holds the ecosystem.

You don't build an ecosphere by giving orders.
You build it by:
  1. Adding rocks (substrate)
  2. Adding soil (nutrients)
  3. Filling with water (event stream)
  4. Planting flora (producers)
  5. Introducing microbiome (decomposers)
  6. Sealing the lid
  7. Adding sunlight (input)

Then you wait for equilibrium.
"""

from __future__ import annotations

import logging
import time
import uuid
from pathlib import Path

from darwin.flora.aragorn import AragornPlant
from darwin.flora.barrnap import BarrnapPlant
from darwin.flora.base import Organism
from darwin.flora.crispard import CRISPARd
from darwin.flora.amr_plant import ABRicatePlant
from darwin.flora.gecco_plant import GeccoPlant
from darwin.flora.isescan_plant import ISEScanPlant
from darwin.flora.minigene import MiniGeneHunter
from darwin.flora.mob_suite_plant import MobSuitePlant
from darwin.flora.phispy_plant import PhiSpyPlant
from darwin.flora.operons import OperonGrouper
from darwin.flora.phylo_16s import PhyloIdentifier
from darwin.flora.prodigal import ProdigalPlant
from darwin.flora.pyhmmer_plant import PyhmmerPlant
from darwin.flora.signal_scanner import SignalScanner
from darwin.microbiome.enricher import Enricher
from darwin.microbiome.scrutinizer import Scrutinizer
from darwin.microbiome.synthesizer import Synthesizer
from darwin.rocks.fasta import parse_fasta
from darwin.rocks.models import AnnotationConfig, Genome
from darwin.soil.nutrients import NutrientStore
from darwin.water.cycle import WaterCycle
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.jar")


class Ecosphere:
    """
    The sealed jar. The complete ecosystem.

    Usage:
        jar = Ecosphere(config)
        result = await jar.add_sunlight(fasta_path)
        # That's it. Equilibrium happens naturally.
    """

    def __init__(self, config: AnnotationConfig | None = None) -> None:
        # The water — everything flows through this
        self.stream = Stream()

        # The soil — databases and tools
        hmm_dbs = config.hmm_databases if config else []
        self.soil = NutrientStore(hmm_databases=hmm_dbs)

        # Configuration
        self.config = config

        # Water cycle observer
        self._cycle: WaterCycle | None = None

        # The organisms — they'll be planted when the jar is sealed
        self._flora: list[Organism] = []
        self._microbiome: list[Organism] = []
        self._sealed = False

        # Fill the jar with life
        self._fill()

    def _fill(self) -> None:
        """
        Fill the jar with organisms.

        Each organism is created and connected to the water
        and soil, but doesn't do anything until sunlight arrives.
        """
        logger.info("🫙 Filling the jar...")

        # Survey the soil first
        soil_report = self.soil.survey()
        logger.info(f"🌱 Soil survey: {soil_report}")

        # Plant the flora (producers)
        self._flora = [
            # Primary producers (feed on GENOME_LOADED)
            ProdigalPlant(self.stream, self.soil),
            AragornPlant(self.stream, self.soil),
            BarrnapPlant(self.stream, self.soil),
            CRISPARd(self.stream, self.soil),
            MobSuitePlant(self.stream, self.soil),
            ISEScanPlant(self.stream, self.soil),
            ABRicatePlant(self.stream, self.soil),
            PhiSpyPlant(self.stream, self.soil),
            GeccoPlant(self.stream, self.soil),
            # Secondary producers (feed on GENES_CALLED)
            PyhmmerPlant(self.stream, self.soil),
            MiniGeneHunter(self.stream, self.soil),
            SignalScanner(self.stream, self.soil),
            OperonGrouper(self.stream, self.soil),
            # Tertiary producers (feed on RRNA_DETECTED)
            PhyloIdentifier(self.stream, self.soil),
        ]

        # Introduce the microbiome (decomposers)
        output_dir = self.config.output_dir if self.config else Path("darwin_output")
        self._microbiome = [
            Scrutinizer(self.stream, self.soil),
            Enricher(self.stream, self.soil),
            Synthesizer(self.stream, self.soil, output_dir=output_dir),
        ]

        # Plant all organisms — they extend roots into the water
        for organism in self._flora + self._microbiome:
            organism.plant()
            logger.info(f"  🌱 Planted {organism.name}")

        # Seal the lid
        self._sealed = True
        logger.info("🫙 Jar sealed. Waiting for sunlight.")

    async def add_sunlight(
        self,
        fasta_path: Path | None = None,
        genome: Genome | None = None,
    ) -> dict:
        """
        Add sunlight to the jar — this starts everything.

        Provide either a FASTA file path or a pre-parsed Genome.
        The sunlight energy flows through the water, waking up
        organisms that feed on it. They produce, decompose,
        and the system reaches equilibrium naturally.

        Returns the final summary when equilibrium is reached.
        """
        if not self._sealed:
            raise RuntimeError("Jar isn't sealed yet! Call _fill() first.")

        # Create a water cycle tracker
        correlation_id = str(uuid.uuid4())[:8]
        self._cycle = WaterCycle(correlation_id=correlation_id)

        # Observe all nutrients flowing
        self.stream.subscribe_all(self._observe)

        # Parse the genome if given a file
        if fasta_path and not genome:
            min_len = self.config.min_contig_length if self.config else 200
            genome = parse_fasta(Path(fasta_path), min_length=min_len)

        if not genome:
            raise ValueError("No sunlight! Provide a FASTA file or Genome.")

        # Build config dict for passing through water
        config_dict = {}
        if self.config:
            config_dict = {
                "metagenome_mode": self.config.metagenome_mode,
                "translation_table": self.config.translation_table,
                "locus_tag_prefix": self.config.locus_tag_prefix,
                "evalue_threshold": self.config.evalue_threshold,
                "output_dir": str(self.config.output_dir),
                "cpus": self.config.cpus,
            }

        logger.info(f"☀️ Sunlight entering the jar: {genome.name}")
        start_time = time.time()

        # Release the genome into the water — this triggers EVERYTHING
        await self.stream.release(
            Nutrient(
                type=NutrientType.GENOME_LOADED,
                data={"genome": genome, "config": config_dict},
                source="sunlight",
                correlation_id=correlation_id,
            )
        )

        # Wait for equilibrium
        reached = await self.stream.wait_for_equilibrium(timeout=600)
        elapsed = time.time() - start_time

        if reached:
            logger.info(f"🌊 Equilibrium reached in {elapsed:.1f}s")
        else:
            logger.warning(f"⚠️ Equilibrium not reached after {elapsed:.1f}s")

        # Collect results from the sediment
        results = self._collect_results()
        results["correlation_id"] = correlation_id
        results["duration_seconds"] = round(elapsed, 2)
        results["equilibrium"] = reached
        results["water_cycle"] = self._cycle.summary() if self._cycle else {}

        return results

    async def _observe(self, nutrient: Nutrient) -> None:
        """Observe nutrients flowing — just watching, not interfering."""
        if self._cycle:
            self._cycle.record(nutrient)

        # Save checkpoints at key stages
        checkpoint_stages = {
            NutrientType.GENES_CALLED: "genes_called",
            NutrientType.PROTEINS_FOUND: "proteins_found",
            NutrientType.QC_COMPLETED: "qc_completed",
        }
        if nutrient.type in checkpoint_stages and "genome" in (nutrient.data or {}):
            output_dir = self.config.output_dir if self.config else Path("darwin_output")
            try:
                from darwin.jar.checkpoint import save_checkpoint

                save_checkpoint(
                    genome=nutrient.data["genome"],
                    stage=checkpoint_stages[nutrient.type],
                    output_dir=output_dir,
                    config=nutrient.data.get("config", {}),
                    metadata={"source": nutrient.source},
                )
            except Exception as e:
                logger.debug(f"Checkpoint save failed (non-critical): {e}")

    def _collect_results(self) -> dict:
        """Read the sediment to find what was produced."""
        results: dict = {
            "output_files": [],
            "genome_summary": {},
            "qc": {},
            "enrichment": {},
        }

        # Find output.written nutrient in sediment
        for nutrient in self.stream.get_sediment(NutrientType.OUTPUT_WRITTEN):
            results["output_files"] = nutrient.data.get("files", [])
            results["genome_summary"] = nutrient.data.get("summary", {}).get("genome", {})

        # Find QC results
        for nutrient in self.stream.get_sediment(NutrientType.QC_COMPLETED):
            results["qc"] = {
                "passed": nutrient.data.get("passed", 0),
                "total": nutrient.data.get("total", 0),
                "healthy": nutrient.data.get("healthy", False),
                "checks": nutrient.data.get("checks", []),
            }

        # Find enrichment
        for nutrient in self.stream.get_sediment(NutrientType.ANNOTATION_READY):
            results["enrichment"] = nutrient.data.get("enrichment", {})

        # Collect any errors
        errors = self.stream.get_sediment(NutrientType.ERROR)
        results["errors"] = [n.data for n in errors]

        # Nutrient flow summary
        results["nutrient_flow"] = self.stream.get_sediment_summary()

        return results

    def get_soil_report(self) -> dict:
        """Check what's in the soil before adding sunlight."""
        return self.soil.survey()

    def reset(self) -> None:
        """Drain and refill — fresh jar."""
        self.stream.reset()
        self._cycle = None
        self._sealed = False
        self._flora.clear()
        self._microbiome.clear()
        self._fill()
