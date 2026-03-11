"""Page Rational — the decision-making agent.

Page Rational handles logical reasoning and decision-making for the council:
  - Decides whether to proceed despite QC warnings
  - Chooses between conflicting annotations
  - Determines optimal parameters based on genome characteristics
  - Applies heuristics for edge cases
  - Can call for a council vote on contentious decisions

Named "Page Rational" because it brings rational, page-by-page analysis
to every decision the council makes.
"""

from __future__ import annotations

from typing import Any

from darwin.agents.base import Agent, AgentMessage, AgentStatus, MessageType


class PageRational(Agent):
    name = "page_rational"
    role = "Decision-making agent — applies rational analysis to council choices"

    def __init__(self) -> None:
        super().__init__()
        self._pending_decisions: list[dict] = []

    def handle_message(self, message: AgentMessage) -> None:
        """Handle decision requests and QC alerts."""
        if message.msg_type == MessageType.ALERT:
            # Evaluate alerts and decide if they're blockers
            severity = message.payload.get("severity", "info")
            if severity == "critical":
                self._pending_decisions.append({
                    "type": "critical_alert",
                    "from": message.sender,
                    "payload": message.payload,
                })

    def validate_input(self, **kwargs: Any) -> bool:
        return self.context.genome is not None

    def execute(self, **kwargs: Any) -> dict[str, Any]:
        self.set_status(AgentStatus.WORKING)
        self.log.info("[bold yellow]Page Rational[/] — analyzing decisions")

        genome = self.context.get_genome()
        decisions: list[dict] = []

        # Decision 1: Should we proceed with incomplete tool set?
        annotator_result = self.context.get_result("annotator")
        if annotator_result:
            skipped = annotator_result.get("tools_skipped", [])
            if skipped:
                decision = self._decide_proceed_with_missing_tools(skipped)
                decisions.append(decision)

        # Decision 2: Parameter optimization based on genome characteristics
        param_decision = self._optimize_parameters(genome)
        decisions.append(param_decision)

        # Decision 3: Evaluate QC results and decide if output is publishable
        scrutinizer_result = self.context.get_result("scrutinizer")
        if scrutinizer_result:
            quality_decision = self._evaluate_quality(scrutinizer_result)
            decisions.append(quality_decision)

        # Record decisions in context
        for d in decisions:
            self.context.record_decision(
                topic=d["topic"],
                votes={self.name: d["vote"]},
                rationale={self.name: d["rationale"]},
                outcome=d["outcome"],
            )

        result = {
            "decisions_made": len(decisions),
            "decisions": decisions,
            "proceed": all(d["outcome"] != "halt" for d in decisions),
        }

        self.context.store_result(self.name, result)
        self.set_status(AgentStatus.DONE)
        return result

    def _decide_proceed_with_missing_tools(self, skipped: list[str]) -> dict:
        """Decide whether to proceed when some tools are unavailable."""
        critical_tools = {"prodigal"}  # Can't annotate without gene calling
        missing_critical = set(skipped) & critical_tools

        if missing_critical:
            return {
                "topic": "proceed_without_critical_tools",
                "vote": "halt",
                "rationale": f"Critical tools missing: {', '.join(missing_critical)}. "
                             f"Cannot produce meaningful annotation without gene calling.",
                "outcome": "halt",
            }

        return {
            "topic": "proceed_with_missing_tools",
            "vote": "proceed",
            "rationale": f"Non-critical tools skipped: {', '.join(skipped)}. "
                         f"Core annotation can still proceed.",
            "outcome": "proceed",
        }

    def _optimize_parameters(self, genome) -> dict:
        """Suggest parameter optimizations based on genome characteristics."""
        suggestions: list[str] = []

        # Small genome might be incomplete → suggest metagenomic mode
        if genome.total_length < 500_000 and not self.context.config.metagenome:
            suggestions.append(
                "Genome <500kb — consider --metagenome mode if this is a MAG"
            )

        # Very high GC → might need adjusted parameters
        if genome.gc_content > 0.70:
            suggestions.append(
                "Very high GC (>70%) — some gene callers may underpredict"
            )

        # Fragmented assembly
        if genome.num_contigs > 200:
            suggestions.append(
                "Highly fragmented — genes split across contigs may be missed"
            )

        return {
            "topic": "parameter_optimization",
            "vote": "proceed",
            "rationale": "; ".join(suggestions) if suggestions else "Parameters look good",
            "outcome": "proceed",
            "suggestions": suggestions,
        }

    def _evaluate_quality(self, qc_result: dict) -> dict:
        """Evaluate overall quality and decide if results are publishable."""
        if not qc_result.get("passed", True):
            failed_checks = [
                c["name"] for c in qc_result.get("checks", [])
                if c.get("severity") == "critical"
            ]
            if failed_checks:
                return {
                    "topic": "output_quality",
                    "vote": "warn",
                    "rationale": f"Critical QC failures: {', '.join(failed_checks)}. "
                                 f"Results should be reviewed before use.",
                    "outcome": "proceed_with_warnings",
                }

        return {
            "topic": "output_quality",
            "vote": "approve",
            "rationale": "QC checks passed — results are suitable for downstream analysis",
            "outcome": "proceed",
        }
