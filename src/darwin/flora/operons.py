"""
OperonGrouper — Operon clustering organism.

Feeds on: genes.called (predicted CDS)
Produces: operons.grouped (genes clustered into operons)

Like a social ecologist mapping communities — groups nearby,
co-directional genes into predicted operons based on
intergenic distance and strand orientation.

Algorithm: sort CDS by contig+strand+start, cluster adjacent
co-directional genes with intergenic gap ≤ 300bp.
"""

from __future__ import annotations

import logging

from darwin.flora.base import Organism
from darwin.rocks.models import FeatureType, Genome
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.flora.operons")

MAX_INTERGENIC_GAP = 300  # bp between adjacent genes in an operon


class OperonGrouper(Organism):
    """Operon predictor — the community mapper."""

    name = "operon_grouper"
    feeds_on_nutrients = [NutrientType.GENES_CALLED]
    produces_nutrients = [NutrientType.OPERONS_GROUPED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return True  # Pure Python

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """Cluster adjacent co-directional genes into operons."""
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})

        self.logger.info("🧬 Grouping genes into operons...")

        operon_count = 0
        genes_in_operons = 0
        operon_id = 0

        for contig in genome.contigs:
            cds = [f for f in contig.features if f.type == FeatureType.CDS]
            if len(cds) < 2:
                continue

            # Sort by start position
            cds.sort(key=lambda f: f.start)

            # Group adjacent co-directional genes
            current_operon = [cds[0]]

            for i in range(1, len(cds)):
                prev = cds[i - 1]
                curr = cds[i]

                # Same strand and close enough?
                gap = curr.start - prev.end
                same_strand = curr.strand == prev.strand

                if same_strand and 0 <= gap <= MAX_INTERGENIC_GAP:
                    current_operon.append(curr)
                else:
                    # Emit operon if ≥ 2 genes
                    if len(current_operon) >= 2:
                        operon_id += 1
                        operon_tag = f"OP_{operon_id:04d}"
                        for pos, gene in enumerate(current_operon, 1):
                            operon_note = (
                                f"operon={operon_tag};"
                                f"operon_pos={pos}/{len(current_operon)}"
                            )
                            if gene.note:
                                gene.note += f";{operon_note}"
                            else:
                                gene.note = operon_note
                        operon_count += 1
                        genes_in_operons += len(current_operon)

                    current_operon = [curr]

            # Don't forget the last group
            if len(current_operon) >= 2:
                operon_id += 1
                operon_tag = f"OP_{operon_id:04d}"
                for pos, gene in enumerate(current_operon, 1):
                    operon_note = (
                        f"operon={operon_tag};"
                        f"operon_pos={pos}/{len(current_operon)}"
                    )
                    if gene.note:
                        gene.note += f";{operon_note}"
                    else:
                        gene.note = operon_note
                operon_count += 1
                genes_in_operons += len(current_operon)

        self.logger.info(
            f"🧬 Found {operon_count} operons containing "
            f"{genes_in_operons} genes"
        )

        return Nutrient(
            type=NutrientType.OPERONS_GROUPED,
            data={
                "genome": genome,
                "operon_count": operon_count,
                "genes_in_operons": genes_in_operons,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )
