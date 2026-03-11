"""Deep Research Bob — database and literature lookup agent.

Bob provides biological context by querying external databases:
  - NCBI taxonomy for organism identification
  - UniProt for protein characterization
  - KEGG for pathway mapping (future)
  - PubMed for relevant literature (future)

Bob is consulted by the Scrutinizer and Page Rational when they
need external context to make decisions.
"""

from __future__ import annotations

from typing import Any

from darwin.agents.base import Agent, AgentMessage, AgentStatus, MessageType


class DeepResearchBob(Agent):
    name = "deep_research_bob"
    role = "Queries external databases for biological context and validation"

    def __init__(self) -> None:
        super().__init__()
        self._research_cache: dict[str, Any] = {}

    def handle_message(self, message: AgentMessage) -> None:
        """Handle research requests from other agents."""
        if message.msg_type == MessageType.REQUEST:
            query_type = message.payload.get("query_type")
            if query_type == "taxonomy":
                result = self._lookup_taxonomy(message.payload)
                self.send(message.sender, MessageType.RESPONSE, result)
            elif query_type == "protein":
                result = self._lookup_protein(message.payload)
                self.send(message.sender, MessageType.RESPONSE, result)

    def validate_input(self, **kwargs: Any) -> bool:
        return self.context.genome is not None

    def execute(self, **kwargs: Any) -> dict[str, Any]:
        self.set_status(AgentStatus.WORKING)
        self.log.info("[bold blue]Deep Research Bob[/] — gathering biological context")

        genome = self.context.get_genome()

        # Gather basic genome context
        context_report = {
            "genome_size_category": self._categorize_genome_size(genome.total_length),
            "gc_interpretation": self._interpret_gc(genome.gc_content),
            "estimated_proteins": sum(
                1 for f in genome.all_features
                if f.feature_type.value == "CDS"
            ),
            "hypothetical_ratio": self._hypothetical_ratio(genome),
            "recommendations": self._generate_recommendations(genome),
        }

        self.context.store_result(self.name, context_report)
        self.context.scratch["research_context"] = context_report

        self.broadcast(
            MessageType.DATA,
            {"message": "Research context compiled", "context": context_report},
        )

        self.log.info(f"  Genome category: {context_report['genome_size_category']}")
        self.log.info(f"  GC interpretation: {context_report['gc_interpretation']}")

        self.set_status(AgentStatus.DONE)
        return context_report

    def _categorize_genome_size(self, total_bp: int) -> str:
        if total_bp < 500_000:
            return "very_small (possible obligate symbiont or incomplete)"
        elif total_bp < 1_500_000:
            return "small (possible reduced genome, e.g. Mycoplasma)"
        elif total_bp < 4_000_000:
            return "medium (typical prokaryote)"
        elif total_bp < 8_000_000:
            return "large (e.g. Streptomyces, Myxobacteria)"
        else:
            return "very_large (unusual for prokaryotes — verify assembly)"

    def _interpret_gc(self, gc: float) -> str:
        if gc < 0.30:
            return "very_low_gc (typical of AT-rich organisms: Mycoplasma, some Firmicutes)"
        elif gc < 0.45:
            return "low_gc (common in Firmicutes, Bacteroidetes)"
        elif gc < 0.55:
            return "medium_gc (typical of many Proteobacteria)"
        elif gc < 0.65:
            return "high_gc (typical of Actinobacteria)"
        else:
            return "very_high_gc (typical of Streptomyces, Deinococcus)"

    def _hypothetical_ratio(self, genome) -> float:
        from darwin.models import FeatureType
        cds = [f for f in genome.all_features if f.feature_type == FeatureType.CDS]
        if not cds:
            return 0.0
        hypo = sum(1 for f in cds if f.product == "hypothetical protein")
        return round(hypo / len(cds), 3)

    def _generate_recommendations(self, genome) -> list[str]:
        """Generate recommendations for further analysis."""
        recs: list[str] = []
        if genome.total_length < 500_000:
            recs.append("Run CheckM2 to assess genome completeness")
        if self._hypothetical_ratio(genome) > 0.5:
            recs.append("High hypothetical ratio — consider InterProScan for deeper annotation")
        if genome.num_contigs > 100:
            recs.append("Fragmented assembly — consider re-assembling with long reads")
        recs.append("Run FastANI against GTDB for taxonomic placement")
        recs.append("Run antiSMASH for biosynthetic gene cluster detection")
        return recs

    def _lookup_taxonomy(self, payload: dict) -> dict:
        """Placeholder for NCBI taxonomy lookup."""
        return {"status": "not_implemented", "note": "NCBI taxonomy lookup coming soon"}

    def _lookup_protein(self, payload: dict) -> dict:
        """Placeholder for UniProt lookup."""
        return {"status": "not_implemented", "note": "UniProt lookup coming soon"}
