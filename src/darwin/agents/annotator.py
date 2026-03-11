"""Annotator Agent — runs all annotation tools on the genome.

Responsibilities:
  - Run Barrnap (rRNA detection) first
  - Run Aragorn (tRNA/tmRNA detection)
  - Run Prodigal (CDS gene calling)
  - Run pyhmmer (protein functional annotation)
  - Assign locus tags in a consistent order
  - Report annotation statistics to the council
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from darwin.agents.base import Agent, AgentStatus, MessageType
from darwin.annotators import (
    AragornAnnotator,
    BarrnapAnnotator,
    ProdigalAnnotator,
    PyhmmerAnnotator,
)
from darwin.models import FeatureType, Genome


class Annotator(Agent):
    name = "annotator"
    role = "Runs all annotation tools (gene calling, RNA, functional) on the genome"

    def validate_input(self, **kwargs: Any) -> bool:
        genome = self.context.genome
        if genome is None or genome.num_contigs == 0:
            self.log.error("No genome available — Processor must run first")
            return False
        return True

    def execute(self, **kwargs: Any) -> dict[str, Any]:
        self.set_status(AgentStatus.WORKING)
        self.log.info("[bold green]Annotator[/] — running annotation tools")

        config = self.context.config
        genome = self.context.get_genome()

        # Track which tools ran and which were skipped
        tools_run: list[str] = []
        tools_skipped: list[str] = []

        # 1. Barrnap — rRNA (run first to avoid gene calls over rRNA)
        barrnap = BarrnapAnnotator(config)
        if barrnap.check_dependencies():
            genome = barrnap.run(genome)
            tools_run.append("barrnap")
        else:
            tools_skipped.append("barrnap")

        # 2. Aragorn — tRNA/tmRNA
        aragorn = AragornAnnotator(config)
        if aragorn.check_dependencies():
            genome = aragorn.run(genome)
            tools_run.append("aragorn")
        else:
            tools_skipped.append("aragorn")

        # 3. Prodigal — CDS
        prodigal = ProdigalAnnotator(config)
        if prodigal.check_dependencies():
            genome = prodigal.run(genome)
            tools_run.append("prodigal")
        else:
            tools_skipped.append("prodigal")
            self.send(
                "scrutinizer",
                MessageType.ALERT,
                {"issue": "Prodigal not found — no gene calling performed",
                 "severity": "critical"},
            )

        # 4. pyhmmer — functional annotation
        hmm_dbs = self._find_hmm_databases()
        if hmm_dbs:
            pyhmmer = PyhmmerAnnotator(config, db_paths=hmm_dbs)
            if pyhmmer.check_dependencies():
                genome = pyhmmer.run(genome)
                tools_run.append("pyhmmer")
            else:
                tools_skipped.append("pyhmmer")
        else:
            tools_skipped.append("pyhmmer (no databases)")

        # Sort features by position
        for contig in genome.contigs:
            contig.features.sort(key=lambda f: (f.start, f.end))

        # Update genome in context
        self.context.set_genome(genome)

        # Compile statistics
        features = genome.all_features
        stats = {
            "tools_run": tools_run,
            "tools_skipped": tools_skipped,
            "cds_count": sum(1 for f in features if f.feature_type == FeatureType.CDS),
            "trna_count": sum(1 for f in features if f.feature_type == FeatureType.TRNA),
            "rrna_count": sum(1 for f in features if f.feature_type == FeatureType.RRNA),
            "tmrna_count": sum(1 for f in features if f.feature_type == FeatureType.TMRNA),
            "total_features": len(features),
            "hypothetical_count": sum(
                1 for f in features
                if f.feature_type == FeatureType.CDS
                and f.product == "hypothetical protein"
            ),
        }

        self.context.store_result(self.name, stats)

        # Notify the council
        self.broadcast(
            MessageType.DATA,
            {"message": f"Annotation complete: {stats['total_features']} features found",
             "stats": stats},
        )

        if tools_skipped:
            self.send(
                "scrutinizer",
                MessageType.ALERT,
                {"issue": f"Tools skipped: {', '.join(tools_skipped)}",
                 "severity": "warning"},
            )

        self.log.info(
            f"  Annotation complete: [green]{stats['cds_count']}[/] CDS, "
            f"[green]{stats['trna_count']}[/] tRNA, "
            f"[green]{stats['rrna_count']}[/] rRNA"
        )

        self.set_status(AgentStatus.DONE)
        return stats

    def _find_hmm_databases(self) -> list[Path]:
        """Look for HMM databases in standard locations."""
        search_dirs = [
            self.context.config.output_dir / "db",
            Path.home() / ".local" / "share" / "darwin" / "db",
            Path("/opt/darwin/db"),
        ]
        db_names = ["Pfam-A.hmm", "TIGRFAMs.hmm"]
        found: list[Path] = []
        for d in search_dirs:
            if not d.exists():
                continue
            for name in db_names:
                path = d / name
                if path.exists():
                    found.append(path)
        return found
