"""
AragornPlant — tRNA/tmRNA detection organism.

Feeds on: genome.loaded
Produces: trna.detected

Like a ground cover plant — small, quiet, but essential
for the ecosystem's completeness. Detects transfer RNAs
that are critical for understanding genome completeness.

Grows in parallel with Prodigal and Barrnap — they all
feed on the same sunlight (genome.loaded) simultaneously.
"""

from __future__ import annotations

import asyncio
import logging
import re
import tempfile
from pathlib import Path

from darwin.flora.base import Organism
from darwin.rocks.models import Feature, FeatureType, Genome, Strand
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.flora.aragorn")


class AragornPlant(Organism):
    """Aragorn tRNA detector — the ground cover."""

    name = "aragorn"
    feeds_on_nutrients = [NutrientType.GENOME_LOADED]
    produces_nutrients = [NutrientType.TRNA_DETECTED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return self.soil.has_aragorn

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """Run Aragorn to find tRNA and tmRNA genes."""
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        locus_prefix = config.get("locus_tag_prefix", "DARWIN")

        if not self.can_grow():
            self.logger.warning("🏜️ Aragorn not in soil")
            return Nutrient(
                type=NutrientType.TRNA_DETECTED,
                data={"genome": genome, "trna_count": 0, "tmrna_count": 0, "config": config},
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        self.logger.info(f"🌾 Searching for tRNAs in {genome.name}...")

        with tempfile.TemporaryDirectory(prefix="darwin_aragorn_") as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / "input.fasta"
            output_file = tmp / "aragorn.txt"

            with open(input_fasta, "w") as fh:
                for contig in genome.contigs:
                    fh.write(f">{contig.id}\n{contig.sequence}\n")

            cmd = [
                "aragorn",
                "-t",
                "-gcbact",
                "-l",
                "-w",
                "-o",
                str(output_file),
                str(input_fasta),
            ]

            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            _, stderr = await proc.communicate()

            if proc.returncode != 0:
                raise RuntimeError(f"Aragorn failed: {stderr.decode()}")

            features = self._parse_output(output_file, locus_prefix)

            # Attach to genome
            trna_count = 0
            tmrna_count = 0
            for feature in features:
                for contig in genome.contigs:
                    if contig.id == feature.contig_id:
                        contig.features.append(feature)
                        if feature.type == FeatureType.TRNA:
                            trna_count += 1
                        elif feature.type == FeatureType.TMRNA:
                            tmrna_count += 1
                        break

        self.logger.info(f"🌾 Found {trna_count} tRNAs, {tmrna_count} tmRNAs")

        return Nutrient(
            type=NutrientType.TRNA_DETECTED,
            data={
                "genome": genome,
                "trna_count": trna_count,
                "tmrna_count": tmrna_count,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    def _parse_output(self, output_path: Path, locus_prefix: str) -> list[Feature]:
        """Parse Aragorn text output format."""
        features = []
        current_contig = ""
        rna_num = 0

        with open(output_path) as fh:
            for line in fh:
                line = line.strip()

                # Contig header
                if line.startswith(">"):
                    current_contig = line[1:].split()[0]
                    continue

                # tRNA line: "1   tRNA-Ala          c[23,95]     ..."
                match = re.match(
                    r"\s*\d+\s+(tRNA-\w+|tmRNA)\s+([c\[]?)(\[?)(\d+),(\d+)\]",
                    line,
                )
                if not match:
                    # Try alternate format
                    match = re.match(
                        r"\s*\d+\s+(tRNA-\w+|tmRNA)\s+(\w*)\[(\d+),(\d+)\]",
                        line,
                    )
                    if match:
                        rna_type_str = match.group(1)
                        start = int(match.group(3))
                        end = int(match.group(4))
                        is_complement = "c" in match.group(2)
                    else:
                        continue
                else:
                    rna_type_str = match.group(1)
                    start = int(match.group(4))
                    end = int(match.group(5))
                    is_complement = "c" in line.split(rna_type_str)[1].split("[")[0]

                rna_num += 1

                if rna_type_str == "tmRNA":
                    ftype = FeatureType.TMRNA
                    product = "transfer-messenger RNA"
                else:
                    ftype = FeatureType.TRNA
                    product = rna_type_str  # e.g., "tRNA-Ala"

                features.append(
                    Feature(
                        type=ftype,
                        start=start,
                        end=end,
                        strand=Strand.REVERSE if is_complement else Strand.FORWARD,
                        contig_id=current_contig,
                        locus_tag=f"{locus_prefix}_t{rna_num:04d}",
                        product=product,
                        inference="COORDINATES:profile:Aragorn",
                    )
                )

        return features
