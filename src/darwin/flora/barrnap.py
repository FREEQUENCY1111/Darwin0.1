"""
BarrnapPlant — rRNA prediction organism.

Feeds on: genome.loaded
Produces: rrna.detected

Like a tall tree in the ecosystem — it identifies the major
structural RNAs (5S, 16S, 23S) that define the organism's
phylogenetic identity.

Grows alongside Aragorn and Prodigal simultaneously.
"""

from __future__ import annotations

import asyncio
import logging
import tempfile
from pathlib import Path
from typing import Optional

from darwin.flora.base import Organism
from darwin.rocks.models import Genome, Feature, FeatureType, Strand
from darwin.water.stream import Nutrient, NutrientType, Stream
from darwin.soil.nutrients import NutrientStore

logger = logging.getLogger("darwin.flora.barrnap")


class BarrnapPlant(Organism):
    """Barrnap rRNA predictor — the tall tree."""

    name = "barrnap"
    feeds_on_nutrients = [NutrientType.GENOME_LOADED]
    produces_nutrients = [NutrientType.RRNA_DETECTED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return self.soil.has_barrnap

    async def grow(self, nutrient: Nutrient) -> Optional[Nutrient]:
        """Run Barrnap to predict rRNA genes."""
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        locus_prefix = config.get("locus_tag_prefix", "DARWIN")

        if not self.can_grow():
            self.logger.warning("🏜️ Barrnap not in soil")
            return Nutrient(
                type=NutrientType.RRNA_DETECTED,
                data={"genome": genome, "rrna_count": 0, "config": config},
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        self.logger.info(f"🍀 Searching for rRNAs in {genome.name}...")

        with tempfile.TemporaryDirectory(prefix="darwin_barrnap_") as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / "input.fasta"

            with open(input_fasta, "w") as fh:
                for contig in genome.contigs:
                    fh.write(f">{contig.id}\n{contig.sequence}\n")

            cmd = ["barrnap", "--kingdom", "bac", str(input_fasta)]

            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, stderr = await proc.communicate()

            if proc.returncode != 0:
                raise RuntimeError(f"Barrnap failed: {stderr.decode()}")

            features = self._parse_gff(stdout.decode(), locus_prefix)

            # Attach to genome
            rrna_count = 0
            rrna_types: dict[str, int] = {}
            for feature in features:
                for contig in genome.contigs:
                    if contig.id == feature.contig_id:
                        contig.features.append(feature)
                        rrna_count += 1
                        rrna_types[feature.product] = rrna_types.get(feature.product, 0) + 1
                        break

        self.logger.info(f"🍀 Found {rrna_count} rRNAs: {rrna_types}")

        return Nutrient(
            type=NutrientType.RRNA_DETECTED,
            data={
                "genome": genome,
                "rrna_count": rrna_count,
                "rrna_types": rrna_types,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    def _parse_gff(self, gff_text: str, locus_prefix: str) -> list[Feature]:
        """Parse Barrnap GFF3 output from stdout."""
        features = []
        rrna_num = 0

        for line in gff_text.strip().split("\n"):
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue

            contig_id = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            score = float(parts[5]) if parts[5] != "." else 0.0
            strand = Strand.FORWARD if parts[6] == "+" else Strand.REVERSE

            # Parse product from attributes
            attrs = parts[8]
            product = "rRNA"
            name_match = None
            for attr in attrs.split(";"):
                if attr.startswith("Name="):
                    product = attr.split("=", 1)[1]
                elif attr.startswith("name="):
                    product = attr.split("=", 1)[1]

            # Standardize product names
            if "5S" in product or "5s" in product:
                product = "5S ribosomal RNA"
            elif "16S" in product or "16s" in product:
                product = "16S ribosomal RNA"
            elif "23S" in product or "23s" in product:
                product = "23S ribosomal RNA"

            rrna_num += 1
            features.append(Feature(
                type=FeatureType.RRNA,
                start=start,
                end=end,
                strand=strand,
                score=score,
                contig_id=contig_id,
                locus_tag=f"{locus_prefix}_r{rrna_num:04d}",
                product=product,
                inference="COORDINATES:profile:Barrnap",
            ))

        return features
