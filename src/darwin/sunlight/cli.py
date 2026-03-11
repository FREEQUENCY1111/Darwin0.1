"""
CLI — The sun. Direct, warm, powerful.

The primary way to interact with Darwin.
"""

from __future__ import annotations

import asyncio
import sys
from pathlib import Path

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from darwin.jar.ecosphere import Ecosphere
from darwin.rocks.models import AnnotationConfig

console = Console()


@click.group()
@click.version_option()
def cli() -> None:
    """Darwin — Prokaryotic Genome Annotator.

    A self-sustaining ecosphere that annotates your genome.
    Just add sunlight.
    """
    pass


@cli.command()
@click.argument("fasta", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default="darwin_output",
    help="Output directory",
)
@click.option("-p", "--prefix", default="DARWIN", help="Locus tag prefix")
@click.option("-t", "--table", default=11, type=int, help="Translation table (default: 11)")
@click.option("-e", "--evalue", default=1e-10, type=float, help="E-value threshold for HMM search")
@click.option("--min-length", default=200, type=int, help="Minimum contig length")
@click.option("--cpus", default=1, type=int, help="Number of CPUs")
@click.option(
    "--hmm-db", multiple=True, type=click.Path(exists=True, path_type=Path), help="HMM database(s)"
)
@click.option("--metagenome", is_flag=True, help="Metagenome mode (Prodigal -p meta)")
def annotate(
    fasta: Path,
    output: Path,
    prefix: str,
    table: int,
    evalue: float,
    min_length: int,
    cpus: int,
    hmm_db: tuple,
    metagenome: bool,
) -> None:
    """Annotate a prokaryotic genome.

    Provide a FASTA file and Darwin will fill the jar,
    add sunlight, and produce annotated output.

    \b
    Example:
        darwin annotate genome.fasta -o results/ -p MYORG
    """
    config = AnnotationConfig(
        input_file=fasta,
        output_dir=output,
        locus_tag_prefix=prefix,
        translation_table=table,
        evalue_threshold=evalue,
        min_contig_length=min_length,
        cpus=cpus,
        hmm_databases=list(hmm_db),
        metagenome_mode=metagenome,
    )

    console.print(
        Panel.fit(
            "[bold green]Darwin Ecosphere[/bold green]\n"
            f"Input: {fasta}\n"
            f"Output: {output}\n"
            f"Prefix: {prefix}",
            title="☀️ Adding Sunlight",
            border_style="green",
        )
    )

    # Build the jar and add sunlight
    jar = Ecosphere(config)

    # Show soil report
    soil = jar.get_soil_report()
    soil_table = Table(title="🌱 Soil Survey")
    soil_table.add_column("Tool", style="cyan")
    soil_table.add_column("Available", style="green")
    for tool, available in soil.items():
        icon = "✅" if available else "❌"
        soil_table.add_row(tool, icon)
    console.print(soil_table)

    # Run the ecosphere
    try:
        result = asyncio.run(jar.add_sunlight(fasta_path=fasta))
    except Exception as e:
        console.print(f"[bold red]☠️ Ecosystem collapsed: {e}[/bold red]")
        sys.exit(1)

    # Display results
    _display_results(result)


@cli.command()
@click.option("--host", default="0.0.0.0", help="Host to bind to")
@click.option("--port", default=8000, type=int, help="Port to listen on")
def serve(host: str, port: int) -> None:
    """Start the Darwin API server (artificial sunlight)."""
    import uvicorn

    from darwin.sunlight.api import create_app

    console.print(
        Panel.fit(
            f"[bold yellow]API Server[/bold yellow]\nhttp://{host}:{port}",
            title="💡 Artificial Sunlight",
            border_style="yellow",
        )
    )

    app = create_app()
    uvicorn.run(app, host=host, port=port)


@cli.command()
def check() -> None:
    """Check what's in the soil (available tools)."""
    from darwin.soil.nutrients import NutrientStore

    soil = NutrientStore()
    report = soil.survey()

    table = Table(title="🌱 Soil Survey — Available Tools")
    table.add_column("Tool", style="cyan")
    table.add_column("Status", style="green")

    for tool, available in report.items():
        status = "[green]✅ Available[/green]" if available else "[red]❌ Not found[/red]"
        table.add_row(tool, status)

    console.print(table)

    if soil.is_fertile:
        console.print("\n[green]✅ Soil is fertile — ready for sunlight![/green]")
    else:
        console.print("\n[red]❌ Soil is barren — install Prodigal at minimum[/red]")


def _display_results(result: dict) -> None:
    """Pretty-print ecosphere results."""
    console.print()

    # Summary
    genome = result.get("genome_summary", {})
    if genome:
        summary_table = Table(title="🪨 Genome Summary")
        summary_table.add_column("Metric", style="cyan")
        summary_table.add_column("Value", style="white")
        for key, val in genome.items():
            if key == "features_by_type":
                val = ", ".join(f"{k}: {v}" for k, v in val.items())
            summary_table.add_row(key, str(val))
        console.print(summary_table)

    # QC
    qc = result.get("qc", {})
    if qc and qc.get("checks"):
        qc_table = Table(title="🔬 Ecosystem Health")
        qc_table.add_column("Check", style="cyan")
        qc_table.add_column("Status")
        qc_table.add_column("Value", style="white")
        qc_table.add_column("Expected", style="dim")

        for check in qc["checks"]:
            icon = "✅" if check["passed"] else "⚠️"
            qc_table.add_row(check["name"], icon, str(check["value"]), check["expected"])
        console.print(qc_table)

    # Enrichment insights
    enrichment = result.get("enrichment", {})
    if enrichment and enrichment.get("insights"):
        console.print(
            Panel(
                "\n".join(f"  • {i}" for i in enrichment["insights"]),
                title="🧬 Enrichment Insights",
                border_style="blue",
            )
        )

    # Output files
    files = result.get("output_files", [])
    if files:
        console.print(
            Panel(
                "\n".join(f"  📄 {f}" for f in files),
                title="📋 Output Files",
                border_style="green",
            )
        )

    # Nutrient flow
    flow = result.get("nutrient_flow", {})
    if flow:
        flow_str = " → ".join(f"{k}" for k in flow.keys())
        console.print(f"\n💧 [dim]Water cycle: {flow_str}[/dim]")

    duration = result.get("duration_seconds", 0)
    eq = "🌊 Equilibrium" if result.get("equilibrium") else "⚠️ Incomplete"
    console.print(f"\n{eq} reached in [bold]{duration}s[/bold]")

    errors = result.get("errors", [])
    if errors:
        console.print(f"\n[red]☠️ {len(errors)} toxins detected in the water[/red]")
        for err in errors:
            console.print(f"  [red]• {err}[/red]")
