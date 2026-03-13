"""
MobSuitePlant — Plasmid detection and typing organism.

Feeds on: genome.loaded (raw sequence data)
Produces: plasmids.classified (contigs enriched with replicon metadata)

Like a geologist classifying rock strata — this organism examines
each contig and determines whether it's chromosomal or plasmid DNA,
then types the plasmid by its mobility and replicon class.

Uses MOB-suite (mob_recon) for plasmid reconstruction and typing.
Install: conda install -c bioconda mob_suite
"""

from __future__ import annotations

import asyncio
import csv
import logging
import tempfile
from pathlib import Path

from darwin.flora.base import Organism
from darwin.rocks.models import Genome
from darwin.soil.nutrients import NutrientStore
from darwin.water import Stream
from darwin.water.stream import Nutrient, NutrientType

logger = logging.getLogger("darwin.flora.mob_suite")


class MobSuitePlant(Organism):
    """MOB-suite plasmid classifier — the geologist."""

    name = "mob_suite"
    feeds_on_nutrients = [NutrientType.GENOME_LOADED]
    produces_nutrients = [NutrientType.PLASMIDS_CLASSIFIED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return self.soil.has_mob_suite

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """
        Run MOB-suite mob_recon on the genome.

        Classifies each contig as chromosome or plasmid, and
        annotates plasmids with replicon type and mobility class.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        cpus = config.get("cpus", 1)

        if not self.can_grow():
            self.logger.warning("🏜️ mob_recon not in soil — cannot classify plasmids")
            return None

        self.logger.info(f"🧬 Classifying replicons in {genome.name}...")

        with tempfile.TemporaryDirectory(prefix="darwin_mobsuite_") as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / "input.fasta"
            output_dir = tmp / "mob_output"

            # Write genome to temp FASTA
            with open(input_fasta, "w") as fh:
                for contig in genome.contigs:
                    fh.write(f">{contig.id}\n{contig.sequence}\n")

            # Run mob_recon
            cmd = [
                "mob_recon",
                "--infile", str(input_fasta),
                "--outdir", str(output_dir),
                "--num_threads", str(cpus),
                "--force",
            ]

            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, stderr = await proc.communicate()

            if proc.returncode != 0:
                err_msg = stderr.decode()[:500]
                self.logger.warning(
                    f"⚠️ mob_recon failed (exit {proc.returncode}): {err_msg}"
                )
                # Non-fatal — return nutrient with unclassified contigs
                return Nutrient(
                    type=NutrientType.PLASMIDS_CLASSIFIED,
                    data={
                        "genome": genome,
                        "plasmid_count": 0,
                        "classification_failed": True,
                        "config": config,
                    },
                    source=self.name,
                    correlation_id=nutrient.correlation_id,
                )

            # Parse contig_report.txt
            plasmid_count = self._parse_contig_report(
                genome, output_dir / "contig_report.txt"
            )

        self.logger.info(
            f"🧬 Classified {plasmid_count} plasmid(s), "
            f"{genome.num_contigs - plasmid_count} chromosomal contig(s)"
        )

        return Nutrient(
            type=NutrientType.PLASMIDS_CLASSIFIED,
            data={
                "genome": genome,
                "plasmid_count": plasmid_count,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    @staticmethod
    def _parse_contig_report(genome: Genome, report_path: Path) -> int:
        """
        Parse MOB-suite contig_report.txt and enrich Contig metadata.

        Returns the number of plasmid contigs found.
        """
        if not report_path.exists():
            logger.warning("MOB-suite contig_report.txt not found")
            return 0

        contig_map = {c.id: c for c in genome.contigs}
        plasmid_count = 0

        with open(report_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                contig_id = row.get("sample_id", "").strip()
                if not contig_id:
                    # Some MOB-suite versions use different column names
                    contig_id = row.get("contig_id", "").strip()
                if contig_id not in contig_map:
                    continue

                contig = contig_map[contig_id]
                molecule = row.get("molecule_type", "").strip().lower()

                if "plasmid" in molecule:
                    contig.replicon_type = "plasmid"
                    plasmid_count += 1
                else:
                    contig.replicon_type = "chromosome"

                # Extract typing info
                rep = row.get("rep_type(s)", "").strip()
                if not rep:
                    rep = row.get("rep_type", "").strip()
                if rep and rep != "-":
                    contig.rep_type = rep

                mob = row.get("relaxase_type(s)", "").strip()
                if not mob:
                    mob = row.get("mob_suite_mob_type", "").strip()
                if mob and mob != "-":
                    contig.mob_type = mob

        # Mark any remaining unclassified contigs as chromosomal
        for contig in genome.contigs:
            if not contig.replicon_type:
                contig.replicon_type = "chromosome"

        return plasmid_count
