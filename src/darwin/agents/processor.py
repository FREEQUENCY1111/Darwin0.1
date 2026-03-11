"""Processor Agent — ingests and preprocesses raw genome data.

Responsibilities:
  - Parse input FASTA files
  - Filter contigs by minimum length
  - Detect and handle compressed inputs (.gz)
  - Calculate basic assembly statistics
  - Normalize sequence data (uppercase, remove ambiguous chars)
  - Flag potential issues (very short contigs, unusual GC, etc.)
"""

from __future__ import annotations

import gzip
import shutil
from pathlib import Path
from typing import Any

from darwin.agents.base import Agent, AgentStatus, MessageType
from darwin.models import Genome
from darwin.utils.fasta import parse_fasta


class Processor(Agent):
    name = "processor"
    role = "Ingests raw FASTA data, preprocesses and validates assembly quality"

    def validate_input(self, **kwargs: Any) -> bool:
        input_path = self.context.config.input_path
        if not input_path.exists():
            self.log.error(f"Input file not found: {input_path}")
            return False
        return True

    def execute(self, **kwargs: Any) -> dict[str, Any]:
        self.set_status(AgentStatus.WORKING)
        self.log.info("[bold green]Processor[/] — ingesting genome data")

        config = self.context.config
        input_path = config.input_path

        # Handle gzipped input
        if str(input_path).endswith(".gz"):
            decompressed = self._decompress(input_path)
            input_path = decompressed

        # Parse FASTA
        genome = parse_fasta(input_path, min_len=config.min_contig_len)

        # Normalize sequences
        genome = self._normalize(genome)

        # Calculate stats and flag issues
        issues = self._quality_check(genome)

        # Store genome in shared context
        self.context.set_genome(genome)

        result = {
            "contigs_loaded": genome.num_contigs,
            "total_bp": genome.total_length,
            "gc_content": round(genome.gc_content, 4),
            "issues": issues,
        }

        self.context.store_result(self.name, result)

        # Notify the council
        self.broadcast(
            MessageType.STATUS,
            {"message": f"Genome loaded: {genome.num_contigs} contigs, {genome.total_length:,} bp"},
        )

        if issues:
            self.send(
                "scrutinizer",
                MessageType.ALERT,
                {"issues": issues, "stage": "preprocessing"},
            )

        self.log.info(
            f"  Loaded [green]{genome.num_contigs}[/] contigs "
            f"({genome.total_length:,} bp, GC={genome.gc_content:.1%})"
        )

        self.set_status(AgentStatus.DONE)
        return result

    def _decompress(self, gz_path: Path) -> Path:
        """Decompress a .gz file to the work dir."""
        out_path = self.context.config.output_dir / ".darwin_tmp" / gz_path.stem
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        self.log.info(f"  Decompressed: {gz_path.name} → {out_path.name}")
        return out_path

    def _normalize(self, genome: Genome) -> Genome:
        """Uppercase all sequences, strip whitespace."""
        for contig in genome.contigs:
            contig.sequence = contig.sequence.upper().strip()
        return genome

    def _quality_check(self, genome: Genome) -> list[str]:
        """Flag potential assembly issues."""
        issues: list[str] = []

        if genome.num_contigs == 0:
            issues.append("CRITICAL: No contigs passed minimum length filter")
            return issues

        if genome.total_length < 100_000:
            issues.append(
                f"WARNING: Very small genome ({genome.total_length:,} bp) — "
                f"may be incomplete"
            )

        if genome.num_contigs > 500:
            issues.append(
                f"WARNING: Highly fragmented assembly ({genome.num_contigs} contigs)"
            )

        gc = genome.gc_content
        if gc < 0.20 or gc > 0.80:
            issues.append(
                f"WARNING: Unusual GC content ({gc:.1%}) — "
                f"verify this is a prokaryotic genome"
            )

        # Check for very short contigs that made it through
        short = sum(1 for c in genome.contigs if c.length < 500)
        if short > 0:
            issues.append(f"INFO: {short} contigs shorter than 500 bp")

        return issues
