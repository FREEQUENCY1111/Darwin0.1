"""Main pipeline runner — orchestrates the annotation workflow.

Order of operations:
  1. Parse input FASTA → Genome
  2. Barrnap (rRNA)     — run first so Prodigal doesn't call genes over rRNA
  3. Aragorn (tRNA/tmRNA)
  4. Prodigal (CDS)
  5. pyhmmer (protein function)
  6. Sort features by position
  7. Write output files
"""

from __future__ import annotations

import time
from pathlib import Path

from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from darwin.annotators import (
    AragornAnnotator,
    BarrnapAnnotator,
    ProdigalAnnotator,
    PyhmmerAnnotator,
)
from darwin.models import AnnotationConfig, Genome
from darwin.output.gff import write_gff3
from darwin.output.genbank import write_genbank
from darwin.output.json_out import write_json
from darwin.utils.fasta import parse_fasta, write_fasta, write_proteins
from darwin.utils.logging import get_logger

console = Console()
log = get_logger("darwin.pipeline")


class DarwinPipeline:
    """Full annotation pipeline for prokaryotic genomes."""

    def __init__(self, config: AnnotationConfig) -> None:
        self.config = config
        self.config.output_dir.mkdir(parents=True, exist_ok=True)

    def run(self) -> Genome:
        """Execute the full pipeline and return the annotated genome."""
        start_time = time.time()

        console.print(
            Panel(
                f"[bold]Darwin Prokaryotic Genome Annotator[/bold]\n"
                f"Input: {self.config.input_path.name}",
                style="blue",
            )
        )

        # 1. Parse input
        log.info("[bold]Parsing input FASTA...[/]")
        genome = parse_fasta(self.config.input_path, self.config.min_contig_len)
        log.info(
            f"  Loaded [green]{genome.num_contigs}[/] contigs "
            f"({genome.total_length:,} bp, GC={genome.gc_content:.1%})"
        )

        if genome.num_contigs == 0:
            log.error("No contigs passed the minimum length filter!")
            return genome

        # 2. Run annotators in order
        annotators = self._build_annotator_chain()
        for annotator in annotators:
            if annotator.check_dependencies():
                genome = annotator.run(genome)
            else:
                log.warning(
                    f"Skipping {annotator.name} — dependencies not found"
                )

        # 3. Sort features by position within each contig
        for contig in genome.contigs:
            contig.features.sort(key=lambda f: (f.start, f.end))

        # 4. Write outputs
        log.info("[bold]Writing output files...[/]")
        self._write_outputs(genome)

        # 5. Summary
        elapsed = time.time() - start_time
        self._print_summary(genome, elapsed)

        return genome

    def _build_annotator_chain(self) -> list:
        """Build the ordered chain of annotators."""
        config = self.config

        chain = [
            BarrnapAnnotator(config),
            AragornAnnotator(config),
            ProdigalAnnotator(config),
        ]

        # pyhmmer needs database paths — look for them in standard locations
        hmm_dbs = self._find_hmm_databases()
        if hmm_dbs:
            chain.append(PyhmmerAnnotator(config, db_paths=hmm_dbs))
        else:
            log.info(
                "  No HMM databases found — skipping functional annotation. "
                "Place Pfam-A.hmm in ~/.local/share/darwin/db/"
            )

        return chain

    def _find_hmm_databases(self) -> list[Path]:
        """Look for HMM databases in standard locations."""
        search_dirs = [
            self.config.output_dir / "db",
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
                    log.info(f"  Found HMM database: {path}")

        return found

    def _write_outputs(self, genome: Genome) -> None:
        """Write all requested output formats."""
        out = self.config.output_dir
        name = genome.name

        formats = self.config.formats
        written: list[str] = []

        if "gff3" in formats:
            path = write_gff3(genome, out / f"{name}.gff3")
            written.append(str(path))

        if "gbk" in formats:
            path = write_genbank(genome, out / f"{name}.gbk")
            written.append(str(path))

        if "json" in formats:
            path = write_json(genome, out / f"{name}.json")
            written.append(str(path))

        if "fna" in formats:
            path = write_fasta(genome.contigs, out / f"{name}.fna")
            written.append(str(path))

        if "faa" in formats:
            # Write protein sequences
            from Bio.Seq import Seq

            proteins = []
            for contig in genome.contigs:
                for f in contig.features:
                    if f.feature_type.value == "CDS":
                        nuc = contig.sequence[f.start - 1 : f.end]
                        if f.strand.value == "-":
                            nuc = str(Seq(nuc).reverse_complement())
                        try:
                            aa = str(
                                Seq(nuc).translate(
                                    table=self.config.translation_table,
                                    to_stop=True,
                                )
                            )
                            proteins.append((f.locus_tag, f.product, aa))
                        except Exception:
                            continue
            if proteins:
                path = write_proteins(proteins, out / f"{name}.faa")
                written.append(str(path))

        for w in written:
            log.info(f"  Wrote: [cyan]{w}[/]")

    def _print_summary(self, genome: Genome, elapsed: float) -> None:
        """Print a nice summary table."""
        summary = genome.summary()

        table = Table(title="Annotation Summary", show_header=False, style="green")
        table.add_column("Metric", style="bold")
        table.add_column("Value", justify="right")

        table.add_row("Genome", summary["name"])
        table.add_row("Total size", f"{summary['total_bp']:,} bp")
        table.add_row("Contigs", str(summary["num_contigs"]))
        table.add_row("GC content", f"{summary['gc_content']:.1%}")
        table.add_row("─" * 20, "─" * 10)
        table.add_row("CDS (proteins)", str(summary["cds_count"]))
        table.add_row("tRNA", str(summary["trna_count"]))
        table.add_row("rRNA", str(summary["rrna_count"]))
        table.add_row("tmRNA", str(summary["tmrna_count"]))
        table.add_row("Total features", str(summary["total_features"]))
        table.add_row("─" * 20, "─" * 10)
        table.add_row("Time", f"{elapsed:.1f}s")

        console.print(table)
