"""The Orchestrator — coordinates the Darwin Agent Council.

The Orchestrator sits at the center of the council (as drawn on the whiteboard).
It:
  1. Initializes all agents and connects them to Context Jimi
  2. Executes the pipeline in the correct order
  3. Handles agent failures gracefully
  4. Consults Page Rational for decisions
  5. Can halt the pipeline if the Scrutinizer vetoes
  6. Produces the final summary

Pipeline execution order:
  Processor → Annotator → Scrutinizer → Page Rational →
  Deep Research Bob → Synthesizer → Code-checker Chris → Uploader
"""

from __future__ import annotations

import time
from typing import Any

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from darwin.agents import (
    Annotator,
    CodeCheckerChris,
    ContextJimi,
    DeepResearchBob,
    PageRational,
    Processor,
    Scrutinizer,
    Synthesizer,
    Uploader,
)
from darwin.agents.base import Agent, AgentStatus
from darwin.models import AnnotationConfig
from darwin.utils.logging import get_logger

console = Console()
log = get_logger("darwin.orchestrator")


class Orchestrator:
    """Coordinates the Darwin Agent Council for genome annotation."""

    def __init__(self, config: AnnotationConfig) -> None:
        self.config = config
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize Context Jimi (shared state)
        self.context = ContextJimi(config)

        # Initialize all agents
        self.processor = Processor()
        self.annotator = Annotator()
        self.scrutinizer = Scrutinizer()
        self.synthesizer = Synthesizer()
        self.uploader = Uploader()
        self.page_rational = PageRational()
        self.deep_research_bob = DeepResearchBob()
        self.code_checker_chris = CodeCheckerChris()

        # The ordered council
        self._pipeline_agents: list[Agent] = [
            self.processor,
            self.annotator,
            self.scrutinizer,
            self.page_rational,
            self.deep_research_bob,
            self.synthesizer,
            self.code_checker_chris,
            self.uploader,
        ]

        # Connect all agents to Context Jimi
        for agent in self._pipeline_agents:
            agent.connect(self.context)

    def run(self) -> dict[str, Any]:
        """Execute the full agent council pipeline."""
        start_time = time.time()

        console.print(Panel(
            "[bold]Darwin Agent Council[/bold]\n"
            f"Input: {self.config.input_path.name}\n"
            f"Agents: {len(self._pipeline_agents)}",
            style="blue",
            title="🧬 Darwin",
        ))

        log.info("[bold]Council assembled — beginning annotation[/]")

        results: dict[str, Any] = {}
        halted = False

        for agent in self._pipeline_agents:
            agent_name = agent.name
            log.info(f"\n{'─' * 50}")

            # Validate inputs
            if not agent.validate_input():
                log.warning(f"  {agent_name}: input validation failed — skipping")
                agent.set_status(AgentStatus.ERROR)
                self.context.store_error(agent_name, "Input validation failed")
                continue

            # Execute the agent
            try:
                result = agent.execute()
                results[agent_name] = result

                # Check if Page Rational says to halt
                if agent_name == "page_rational" and not result.get("proceed", True):
                    log.error(
                        "[red bold]Page Rational halted the pipeline[/] — "
                        "critical issues detected"
                    )
                    halted = True
                    break

            except Exception as exc:
                log.error(f"  {agent_name} failed: {exc}")
                agent.set_status(AgentStatus.ERROR)
                self.context.store_error(agent_name, str(exc))

                # Non-critical agents can fail without halting
                critical_agents = {"processor", "annotator", "synthesizer"}
                if agent_name in critical_agents:
                    log.error(f"  Critical agent {agent_name} failed — halting pipeline")
                    halted = True
                    break

        elapsed = time.time() - start_time

        # Print summary
        self._print_council_summary(results, elapsed, halted)

        return {
            "halted": halted,
            "elapsed_seconds": round(elapsed, 1),
            "agent_results": results,
            "council_status": self.context.get_council_status(),
            "decisions": [
                {"topic": d.topic, "outcome": d.outcome}
                for d in self.context.get_decisions()
            ],
            "genome_summary": (
                self.context.genome.summary()
                if self.context.genome
                else None
            ),
        }

    def _print_council_summary(
        self,
        results: dict,
        elapsed: float,
        halted: bool,
    ) -> None:
        """Print a summary of the council's work."""
        console.print(f"\n{'═' * 50}")

        # Agent status table
        status_table = Table(title="Agent Council Status", show_header=True)
        status_table.add_column("Agent", style="bold")
        status_table.add_column("Role", style="dim")
        status_table.add_column("Status", justify="center")

        status_map = self.context.get_council_status()
        for agent in self._pipeline_agents:
            status = status_map.get(agent.name, "unknown")
            if status == "done":
                style = "[green]DONE[/]"
            elif status == "error":
                style = "[red]ERROR[/]"
            elif status == "working":
                style = "[yellow]WORKING[/]"
            else:
                style = f"[dim]{status}[/]"
            status_table.add_row(agent.name, agent.role[:50], style)

        console.print(status_table)

        # Genome summary if available
        if self.context.genome:
            summary = self.context.genome.summary()
            genome_table = Table(title="Annotation Results", show_header=False)
            genome_table.add_column("Metric", style="bold")
            genome_table.add_column("Value", justify="right")

            genome_table.add_row("Genome", summary["name"])
            genome_table.add_row("Size", f"{summary['total_bp']:,} bp")
            genome_table.add_row("Contigs", str(summary["num_contigs"]))
            genome_table.add_row("GC", f"{summary['gc_content']:.1%}")
            genome_table.add_row("CDS", str(summary["cds_count"]))
            genome_table.add_row("tRNA", str(summary["trna_count"]))
            genome_table.add_row("rRNA", str(summary["rrna_count"]))
            genome_table.add_row("Total features", str(summary["total_features"]))
            genome_table.add_row("Time", f"{elapsed:.1f}s")

            console.print(genome_table)

        # Decisions
        decisions = self.context.get_decisions()
        if decisions:
            console.print(f"\n[bold]Council Decisions:[/]")
            for d in decisions:
                console.print(f"  • {d.topic}: [cyan]{d.outcome}[/]")

        if halted:
            console.print(
                Panel("[red bold]Pipeline HALTED — review issues above[/]", style="red")
            )
        else:
            console.print(
                Panel(
                    f"[green bold]Annotation complete![/]\n"
                    f"Results: {self.config.output_dir}",
                    style="green",
                )
            )
