"""
Enricher — Contextual analysis decomposer.

Feeds on: proteins.found + qc.completed
Produces: annotation.ready

Like nitrogen-fixing bacteria — takes raw material and
enriches it with contextual information. Analyzes the genome
characteristics, hypothetical protein ratio, and generates
research insights.

This is where Deep Research Bob's wisdom lives now —
not as a boss, but as a humble microbe doing its thing.
"""

from __future__ import annotations

import logging
from typing import Optional

from darwin.flora.base import Organism
from darwin.rocks.models import Genome, FeatureType
from darwin.water.stream import Nutrient, NutrientType, Stream
from darwin.soil.nutrients import NutrientStore

logger = logging.getLogger("darwin.microbiome.enricher")


class Enricher(Organism):
    """
    Context enrichment microbe.

    Waits for proteins to be annotated AND QC to pass,
    then enriches the annotation with contextual analysis.
    """

    name = "enricher"
    feeds_on_nutrients = [NutrientType.QC_COMPLETED]
    produces_nutrients = [NutrientType.ANNOTATION_READY]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    async def grow(self, nutrient: Nutrient) -> Optional[Nutrient]:
        """Enrich the annotation with contextual analysis."""
        genome: Genome = nutrient.data["genome"]
        qc_data = nutrient.data
        config = nutrient.data.get("config", {})

        self.logger.info(f"🧬 Enriching annotation for {genome.name}...")

        # Genome characterization
        size_category = self._categorize_size(genome.total_length)
        gc_interpretation = self._interpret_gc(genome.gc_content)

        # Feature statistics
        cds = genome.features_of_type(FeatureType.CDS)
        trnas = genome.features_of_type(FeatureType.TRNA)
        rrnas = genome.features_of_type(FeatureType.RRNA)

        hypothetical = sum(1 for f in cds if f.product == "hypothetical protein")
        hyp_ratio = hypothetical / len(cds) if cds else 0

        # Generate insights
        insights = []

        insights.append(f"Genome classification: {size_category}")
        insights.append(f"GC content analysis: {gc_interpretation}")
        insights.append(
            f"Hypothetical protein ratio: {hyp_ratio:.1%} "
            f"({hypothetical}/{len(cds)})"
        )

        if hyp_ratio > 0.5:
            insights.append(
                "HIGH hypothetical ratio suggests novel organism or "
                "incomplete reference databases"
            )
        elif hyp_ratio > 0.3:
            insights.append(
                "Moderate hypothetical ratio — consider searching "
                "additional databases (UniRef, COG)"
            )

        if len(trnas) < 30:
            insights.append(
                f"Low tRNA count ({len(trnas)}) — genome may be "
                "incomplete or highly reduced"
            )

        # Recommendations
        recommendations = self._generate_recommendations(
            genome, cds, trnas, rrnas, hyp_ratio, qc_data
        )

        self.logger.info(f"🧬 Generated {len(insights)} insights, {len(recommendations)} recommendations")

        return Nutrient(
            type=NutrientType.ANNOTATION_READY,
            data={
                "genome": genome,
                "enrichment": {
                    "size_category": size_category,
                    "gc_interpretation": gc_interpretation,
                    "hypothetical_ratio": round(hyp_ratio, 4),
                    "insights": insights,
                    "recommendations": recommendations,
                    "stats": {
                        "total_cds": len(cds),
                        "total_trna": len(trnas),
                        "total_rrna": len(rrnas),
                        "hypothetical": hypothetical,
                        "annotated": len(cds) - hypothetical,
                    },
                },
                "qc": qc_data.get("checks", []),
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    def _categorize_size(self, total_bp: int) -> str:
        if total_bp < 600_000:
            return "minimal/endosymbiont (<0.6 Mb)"
        elif total_bp < 2_000_000:
            return "small prokaryote (0.6-2 Mb)"
        elif total_bp < 5_000_000:
            return "typical prokaryote (2-5 Mb)"
        elif total_bp < 8_000_000:
            return "large prokaryote (5-8 Mb)"
        else:
            return "very large prokaryote (>8 Mb)"

    def _interpret_gc(self, gc: float) -> str:
        if gc < 30:
            return f"{gc}% — very AT-rich (obligate intracellular/endosymbiont)"
        elif gc < 40:
            return f"{gc}% — AT-rich (Firmicutes-like)"
        elif gc < 50:
            return f"{gc}% — moderate"
        elif gc < 60:
            return f"{gc}% — GC-rich (Actinobacteria-like)"
        elif gc < 70:
            return f"{gc}% — high GC (Streptomyces-like)"
        else:
            return f"{gc}% — extremely GC-rich"

    def _generate_recommendations(
        self,
        genome: Genome,
        cds: list,
        trnas: list,
        rrnas: list,
        hyp_ratio: float,
        qc_data: dict,
    ) -> list[str]:
        recs = []

        if hyp_ratio > 0.4:
            recs.append("Run against UniRef90/UniRef50 for additional annotations")
            recs.append("Consider InterProScan for domain-level analysis")

        if not rrnas:
            recs.append("No rRNAs found — check if Barrnap is installed")
        elif len(rrnas) < 3:
            recs.append("Incomplete rRNA operon — genome may be fragmented")

        if genome.num_contigs > 100:
            recs.append(
                f"Highly fragmented ({genome.num_contigs} contigs) — "
                "consider reassembly or scaffolding"
            )

        if genome.total_length > 8_000_000:
            recs.append("Very large genome — consider antiSMASH for secondary metabolites")

        if qc_data.get("critical_failures", 0) > 0:
            recs.append("Critical QC failures detected — review results carefully")

        if not recs:
            recs.append("Annotation looks good — no additional analysis recommended")

        return recs
