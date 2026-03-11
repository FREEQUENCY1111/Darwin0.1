"""Darwin CLI — annotate prokaryotic genomes from the command line.

Usage:
    darwin annotate genome.fasta -o output/
    darwin annotate genome.fasta --metagenome --cpus 8
    darwin serve --port 8000
"""

from __future__ import annotations

from pathlib import Path

import click
from rich.console import Console

from darwin import __version__

console = Console()


@click.group()
@click.version_option(version=__version__, prog_name="darwin")
def main() -> None:
    """Darwin — Fast, accurate prokaryotic genome annotation."""
    pass


@main.command()
@click.argument("input_fasta", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o", "--output",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory [default: <input_name>_darwin]",
)
@click.option(
    "--locus-tag",
    default="DARWIN",
    show_default=True,
    help="Locus tag prefix for features",
)
@click.option(
    "--translation-table",
    type=int,
    default=11,
    show_default=True,
    help="NCBI translation table (11=bacterial, 4=mycoplasma)",
)
@click.option(
    "--evalue",
    type=float,
    default=1e-6,
    show_default=True,
    help="E-value threshold for HMM searches",
)
@click.option(
    "--cpus",
    type=int,
    default=1,
    show_default=True,
    help="Number of CPUs to use",
)
@click.option(
    "--metagenome",
    is_flag=True,
    default=False,
    help="Use metagenomic mode for gene prediction",
)
@click.option(
    "--kingdom",
    type=click.Choice(["bac", "arc"]),
    default="bac",
    show_default=True,
    help="Kingdom for rRNA prediction",
)
@click.option(
    "--min-contig-len",
    type=int,
    default=200,
    show_default=True,
    help="Minimum contig length to annotate",
)
@click.option(
    "--formats",
    default="gff3,gbk,faa,fna,json",
    show_default=True,
    help="Comma-separated output formats",
)
def annotate(
    input_fasta: Path,
    output: Path | None,
    locus_tag: str,
    translation_table: int,
    evalue: float,
    cpus: int,
    metagenome: bool,
    kingdom: str,
    min_contig_len: int,
    formats: str,
) -> None:
    """Annotate a prokaryotic genome from a FASTA file."""
    from darwin.models import AnnotationConfig
    from darwin.council.orchestrator import Orchestrator

    if output is None:
        output = input_fasta.parent / f"{input_fasta.stem}_darwin"

    config = AnnotationConfig(
        input_path=input_fasta,
        output_dir=output,
        locus_tag_prefix=locus_tag,
        translation_table=translation_table,
        evalue=evalue,
        cpus=cpus,
        metagenome=metagenome,
        kingdom=kingdom,
        min_contig_len=min_contig_len,
        formats=[f.strip() for f in formats.split(",")],
    )

    # Run the Agent Council
    orchestrator = Orchestrator(config)
    result = orchestrator.run()

    if result.get("halted"):
        console.print("\n[red bold]Pipeline halted — see issues above[/]")
        raise SystemExit(1)
    else:
        console.print(f"\n[green bold]Done![/] Results in: {output}")


@main.command()
@click.option("--host", default="0.0.0.0", show_default=True, help="Host to bind")
@click.option("--port", type=int, default=8000, show_default=True, help="Port to bind")
@click.option("--reload", is_flag=True, help="Auto-reload on code changes (dev mode)")
def serve(host: str, port: int, reload: bool) -> None:
    """Start the Darwin REST API server."""
    import uvicorn

    console.print(f"[bold blue]Starting Darwin API on {host}:{port}[/]")
    uvicorn.run(
        "darwin.api.app:app",
        host=host,
        port=port,
        reload=reload,
    )


@main.command()
def check() -> None:
    """Check that all external dependencies are installed."""
    from darwin.annotators import (
        AragornAnnotator,
        BarrnapAnnotator,
        ProdigalAnnotator,
        PyhmmerAnnotator,
    )
    from darwin.models import AnnotationConfig

    # Dummy config just to instantiate annotators
    dummy = AnnotationConfig(
        input_path=Path("/dev/null"),
        output_dir=Path("/tmp/darwin_check"),
    )

    tools = [
        ("Prodigal", ProdigalAnnotator(dummy)),
        ("Barrnap", BarrnapAnnotator(dummy)),
        ("Aragorn", AragornAnnotator(dummy)),
        ("pyhmmer", PyhmmerAnnotator(dummy)),
    ]

    console.print("[bold]Checking Darwin dependencies...[/]\n")
    all_ok = True
    for name, annotator in tools:
        ok = annotator.check_dependencies()
        status = "[green]OK[/]" if ok else "[red]MISSING[/]"
        console.print(f"  {name:15s} {status}")
        if not ok:
            all_ok = False

    if all_ok:
        console.print("\n[green bold]All dependencies available![/]")
    else:
        console.print(
            "\n[yellow]Some tools are missing. Install with:[/]\n"
            "  conda install -c bioconda prodigal barrnap aragorn\n"
            "  pip install pyhmmer"
        )
