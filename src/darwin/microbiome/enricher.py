"""
Enricher — Contextual analysis decomposer.

Feeds on: qc.completed
Produces: annotation.ready

Like nitrogen-fixing bacteria — takes raw material and
enriches it with contextual information. Analyzes the genome
characteristics, hypothetical protein ratio, metabolic
markers, potential HGT regions, and generates research insights.
"""

from __future__ import annotations

import logging

from darwin.flora.base import Organism
from darwin.rocks.models import FeatureType, Genome
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.microbiome.enricher")

# Metabolic marker keywords for lifestyle inference
METABOLIC_MARKERS = {
    "aerobic": [
        "cytochrome c oxidase", "ATP synthase", "NADH dehydrogenase",
        "succinate dehydrogenase", "cytochrome bc1",
    ],
    "anaerobic": [
        "fumarate reductase", "nitrate reductase", "formate dehydrogenase",
        "hydrogenase", "ferredoxin",
    ],
    "nitrogen_fixing": [
        "nitrogenase", "nifH", "nifD", "nifK", "nitrogen fixation",
    ],
    "phototropic": [
        "photosystem", "bacteriochlorophyll", "light-harvesting",
        "reaction center", "RuBisCO",
    ],
    "motile": [
        "flagellin", "flagellar", "chemotaxis", "methyl-accepting",
        "CheA", "CheY", "MotA", "MotB",
    ],
    "sporulating": [
        "sporulation", "spore coat", "SpoIIE", "SpoVA",
        "dipicolinate synthase",
    ],
    "antibiotic_resistance": [
        "beta-lactamase", "aminoglycoside", "tetracycline resistance",
        "chloramphenicol acetyltransferase", "vancomycin resistance",
        "efflux pump", "multidrug resistance",
    ],
    "pathogenicity": [
        "hemolysin", "adhesin", "invasin", "toxin", "type III secretion",
        "type IV secretion", "virulence",
    ],
}


class Enricher(Organism):
    """
    Context enrichment microbe.

    Waits for QC to complete, then enriches the annotation
    with contextual analysis including metabolic markers,
    lifestyle inference, and genomic island detection.
    """

    name = "enricher"
    feeds_on_nutrients = [NutrientType.QC_COMPLETED]
    produces_nutrients = [NutrientType.ANNOTATION_READY]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
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
        crisprs = genome.features_of_type(FeatureType.CRISPR)
        signal_peptides = genome.features_of_type(FeatureType.SIGNAL_PEPTIDE)

        hypothetical = sum(1 for f in cds if f.product == "hypothetical protein")
        hyp_ratio = hypothetical / len(cds) if cds else 0

        # Generate insights
        insights = []

        insights.append(f"Genome classification: {size_category}")
        insights.append(f"GC content analysis: {gc_interpretation}")
        insights.append(f"Hypothetical protein ratio: {hyp_ratio:.1%} ({hypothetical}/{len(cds)})")

        if hyp_ratio > 0.5:
            insights.append(
                "HIGH hypothetical ratio suggests novel organism or incomplete reference databases"
            )
        elif hyp_ratio > 0.3:
            insights.append(
                "Moderate hypothetical ratio — consider searching "
                "additional databases (UniRef, COG)"
            )

        if len(trnas) < 30:
            insights.append(
                f"Low tRNA count ({len(trnas)}) — genome may be incomplete or highly reduced"
            )

        # CRISPR insights
        if crisprs:
            total_spacers = sum(
                int(s.split("spacers")[0].split(",")[-1].strip())
                for s in [f.product for f in crisprs]
                if "spacers" in s
            )
            insights.append(
                f"CRISPR defense: {len(crisprs)} array(s) with "
                f"~{total_spacers} spacers — active adaptive immunity"
            )

        # Signal peptide insights
        if signal_peptides:
            sp_ratio = len(signal_peptides) / len(cds) * 100 if cds else 0
            insights.append(
                f"Secretome: {len(signal_peptides)} signal peptides "
                f"({sp_ratio:.1f}% of proteins)"
            )

        # Operon insights
        operon_genes = sum(1 for f in cds if f.note and "operon=" in f.note)
        if operon_genes:
            operon_ids = set()
            for f in cds:
                if f.note and "operon=" in f.note:
                    for part in f.note.split(";"):
                        if part.startswith("operon="):
                            operon_ids.add(part.split("=")[1])
            insights.append(
                f"Gene organization: {len(operon_ids)} predicted operons "
                f"({operon_genes} genes)"
            )

        # Replicon structure insights
        plasmid_contigs = [c for c in genome.contigs if c.replicon_type == "plasmid"]
        chromo_contigs = [c for c in genome.contigs if c.replicon_type == "chromosome"]
        if plasmid_contigs or chromo_contigs:
            insights.append(
                f"Replicon structure: {len(chromo_contigs)} chromosome(s), "
                f"{len(plasmid_contigs)} plasmid(s)"
            )
            for pc in plasmid_contigs:
                details = [f"{pc.length:,} bp"]
                if pc.rep_type:
                    details.append(f"replicon: {pc.rep_type}")
                if pc.mob_type:
                    details.append(f"mobility: {pc.mob_type}")
                insights.append(f"  Plasmid {pc.id}: {', '.join(details)}")
            # Conjugation potential
            conjugative = [c for c in plasmid_contigs if c.mob_type]
            if conjugative:
                insights.append(
                    f"Conjugation potential: {len(conjugative)} plasmid(s) "
                    f"with mobilization genes"
                )

        # Mobile element insights
        mobile_elements = genome.features_of_type(FeatureType.MOBILE_ELEMENT)
        if mobile_elements:
            families = set()
            for f in mobile_elements:
                if "IS family:" in f.note:
                    fam = f.note.split("IS family:")[1].split(";")[0].strip()
                    families.add(fam)
            insights.append(
                f"Mobile elements: {len(mobile_elements)} IS element(s) "
                f"from {len(families)} family/families"
                + (f" ({', '.join(sorted(families)[:5])})" if families else "")
            )

        # AMR insights
        amr_genes = genome.features_of_type(FeatureType.AMR_GENE)
        if amr_genes:
            resistance_classes: set[str] = set()
            for f in amr_genes:
                if "resistance class:" in f.note:
                    cls_str = f.note.split("resistance class:")[1].split(";")[0].strip()
                    for cls in cls_str.split(";"):
                        resistance_classes.add(cls.strip())
            insights.append(
                f"Antimicrobial resistance: {len(amr_genes)} resistance gene(s)"
                + (f"; drug classes: {', '.join(sorted(resistance_classes)[:6])}"
                   if resistance_classes else "")
            )

        # Prophage insights
        prophages = genome.features_of_type(FeatureType.PROPHAGE)
        if prophages:
            intact = sum(1 for f in prophages if "intact" in f.note.lower())
            insights.append(
                f"Viral integration: {len(prophages)} prophage region(s)"
                + (f"; {intact} intact, {len(prophages) - intact} incomplete/predicted"
                   if intact else "")
            )

        # BGC insights
        bgcs = genome.features_of_type(FeatureType.BGC)
        if bgcs:
            bgc_types: dict[str, int] = {}
            for f in bgcs:
                if "type:" in f.note:
                    btype = f.note.split("type:")[1].split(";")[0].strip()
                    bgc_types[btype] = bgc_types.get(btype, 0) + 1
            type_str = ", ".join(f"{k}({v})" for k, v in sorted(bgc_types.items()))
            insights.append(
                f"Biosynthetic potential: {len(bgcs)} gene cluster(s)"
                + (f"; types: {type_str}" if type_str else "")
            )

        # Functional annotation insights (InterProScan)
        go_proteins = [f for f in cds if f.go_terms]
        ipr_proteins = [f for f in cds if f.ipr_ids]
        if go_proteins or ipr_proteins:
            all_go_terms = set()
            for f in go_proteins:
                all_go_terms.update(f.go_terms)
            all_ipr_ids = set()
            for f in ipr_proteins:
                all_ipr_ids.update(f.ipr_ids)
            insights.append(
                f"Functional annotation: {len(go_proteins)} proteins with GO terms "
                f"({len(all_go_terms)} unique), {len(ipr_proteins)} with InterPro "
                f"domains ({len(all_ipr_ids)} unique)"
            )

        # Structural homology insights (Foldseek)
        struct_proteins = [f for f in cds if f.structure_hit]
        if struct_proteins:
            remote_homologs = [
                f for f in struct_proteins
                if "seq identity" in f.note
                and any(
                    float(p.split(":")[1].strip().rstrip("%")) < 30.0
                    for p in f.note.split(";")
                    if p.strip().startswith("seq identity")
                )
            ]
            insights.append(
                f"Structural homology: {len(struct_proteins)} proteins with "
                f"PDB/AlphaFold hits"
                + (f", {len(remote_homologs)} remote homologs (<30% identity)"
                   if remote_homologs else "")
            )

        # Taxonomy insight
        if genome.taxonomy:
            insights.append(f"Taxonomic classification: {genome.taxonomy}")

        # Metabolic markers
        metabolic = self._scan_metabolic_markers(cds)
        if metabolic:
            lifestyle = self._infer_lifestyle(metabolic)
            insights.append(f"Metabolic capabilities: {', '.join(metabolic.keys())}")
            if lifestyle:
                insights.append(f"Predicted lifestyle: {lifestyle}")

        # GC island detection (putative HGT)
        gc_islands = self._detect_gc_islands(genome)
        if gc_islands:
            insights.append(
                f"Putative genomic islands: {len(gc_islands)} regions with "
                f"aberrant GC content (potential HGT)"
            )

        # Recommendations
        recommendations = self._generate_recommendations(
            genome, cds, trnas, rrnas, hyp_ratio, qc_data
        )

        self.logger.info(
            f"🧬 Generated {len(insights)} insights, {len(recommendations)} recommendations"
        )

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
                    "metabolic_markers": {k: len(v) for k, v in metabolic.items()} if metabolic else {},
                    "gc_islands": len(gc_islands),
                    "stats": {
                        "total_cds": len(cds),
                        "total_trna": len(trnas),
                        "total_rrna": len(rrnas),
                        "total_crispr": len(crisprs),
                        "total_signal_peptides": len(signal_peptides),
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
        elif gc < 55:
            return f"{gc}% — moderate (typical Proteobacteria/Enterobacteriaceae)"
        elif gc < 65:
            return f"{gc}% — GC-rich (Actinobacteria-like)"
        elif gc < 72:
            return f"{gc}% — high GC (Streptomyces-like)"
        else:
            return f"{gc}% — extremely GC-rich"

    def _scan_metabolic_markers(self, cds: list) -> dict[str, list[str]]:
        """Scan product names for metabolic pathway markers."""
        found: dict[str, list[str]] = {}
        for f in cds:
            product_lower = f.product.lower()
            for category, keywords in METABOLIC_MARKERS.items():
                for keyword in keywords:
                    if keyword.lower() in product_lower:
                        if category not in found:
                            found[category] = []
                        found[category].append(f.product)
                        break  # one hit per feature per category
        return found

    def _infer_lifestyle(self, metabolic: dict[str, list[str]]) -> str:
        """Infer organism lifestyle from metabolic markers."""
        traits = []
        if "aerobic" in metabolic:
            traits.append("aerobic respiration")
        if "anaerobic" in metabolic:
            traits.append("anaerobic metabolism")
        if "nitrogen_fixing" in metabolic:
            traits.append("nitrogen-fixing")
        if "phototropic" in metabolic:
            traits.append("phototrophic")
        if "motile" in metabolic:
            traits.append("motile")
        if "sporulating" in metabolic:
            traits.append("spore-forming")
        if "pathogenicity" in metabolic:
            traits.append("potential pathogen")
        if "antibiotic_resistance" in metabolic:
            count = len(metabolic["antibiotic_resistance"])
            traits.append(f"antibiotic resistance ({count} genes)")
        return ", ".join(traits) if traits else ""

    def _detect_gc_islands(self, genome: Genome, window_size: int = 5000) -> list[dict]:
        """
        Detect regions where GC deviates >8% from genome mean.

        These putative genomic islands may represent HGT events.
        """
        islands: list[dict] = []
        genome_gc = genome.gc_content

        for contig in genome.contigs:
            seq = contig.sequence.upper()
            if len(seq) < window_size * 2:
                continue

            for i in range(0, len(seq) - window_size + 1, window_size // 2):
                window = seq[i : i + window_size]
                gc = sum(1 for b in window if b in "GC") / len(window) * 100
                deviation = abs(gc - genome_gc)

                if deviation > 8.0:
                    islands.append({
                        "contig": contig.id,
                        "start": i + 1,
                        "end": i + window_size,
                        "gc": round(gc, 1),
                        "deviation": round(deviation, 1),
                    })

        return islands

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
