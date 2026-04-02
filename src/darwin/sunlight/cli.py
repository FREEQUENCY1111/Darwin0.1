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
@click.option("--resume", is_flag=True, help="Resume from checkpoint if available")
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
    resume: bool,
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

    # Check for resume
    if resume:
        from darwin.jar.checkpoint import list_checkpoints

        checkpoints = list_checkpoints(output)
        if checkpoints:
            latest = checkpoints[-1]
            console.print(
                f"[yellow]💾 Found checkpoint: {latest['stage']} "
                f"({latest['saved_at']})[/yellow]"
            )
            console.print("[yellow]   Resuming from checkpoint...[/yellow]")
        else:
            console.print("[dim]No checkpoints found — starting fresh.[/dim]")

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
@click.option("--pfam", is_flag=True, help="Download Pfam-A database")
@click.option("--tigrfams", is_flag=True, help="Download TIGRFAMs database")
@click.option("--all", "download_all_dbs", is_flag=True, help="Download all databases")
@click.option("--force", is_flag=True, help="Re-download even if cached")
def setup(pfam: bool, tigrfams: bool, download_all_dbs: bool, force: bool) -> None:
    """Download and prepare HMM databases.

    Fetches databases to ~/.darwin/databases/ so they're
    automatically available for annotation runs.

    \b
    Example:
        darwin setup --all
        darwin setup --pfam --tigrfams
    """
    from darwin.soil.downloader import download_all, download_database, get_cache_dir

    console.print(
        Panel.fit(
            f"[bold cyan]Database Setup[/bold cyan]\nCache: {get_cache_dir()}",
            title="🌱 Enriching the Soil",
            border_style="cyan",
        )
    )

    if download_all_dbs or (not pfam and not tigrfams):
        # Default: download all if nothing specified
        if not download_all_dbs and not pfam and not tigrfams:
            console.print("[dim]No database specified — downloading all.[/dim]")
        results = download_all(force=force)
        if results:
            console.print(f"\n[green]✅ {len(results)} database(s) ready![/green]")
        else:
            console.print("\n[red]❌ No databases downloaded.[/red]")
        return

    downloaded = []
    if pfam:
        path = download_database("pfam", force=force)
        if path:
            downloaded.append(path)
            console.print(f"[green]✅ Pfam-A ready: {path}[/green]")
    if tigrfams:
        path = download_database("tigrfams", force=force)
        if path:
            downloaded.append(path)
            console.print(f"[green]✅ TIGRFAMs ready: {path}[/green]")

    if downloaded:
        console.print(f"\n[green]✅ {len(downloaded)} database(s) ready![/green]")
    else:
        console.print("\n[red]❌ No databases downloaded.[/red]")


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


@cli.command()
@click.argument("genome1", type=click.Path(exists=True, path_type=Path))
@click.argument("genome2", type=click.Path(exists=True, path_type=Path))
@click.option("-k", "--kmer-size", default=16, type=int, help="K-mer size (default: 16)")
@click.option(
    "-s", "--sketch-size", default=10000, type=int, help="Sketch size (default: 10000)"
)
def compare(genome1: Path, genome2: Path, kmer_size: int, sketch_size: int) -> None:
    """Compare two genomes and estimate ANI.

    Uses k-mer Jaccard similarity to approximate Average
    Nucleotide Identity. Species threshold: ANI ≥ 95%.

    \b
    Example:
        darwin compare genome1.fna genome2.fna
    """
    from darwin.utils.ani import compare_genomes

    console.print(
        Panel.fit(
            f"[bold cyan]Genome Comparison[/bold cyan]\n"
            f"Genome 1: {genome1}\n"
            f"Genome 2: {genome2}\n"
            f"k={kmer_size}, sketch={sketch_size}",
            title="🔬 ANI Estimation",
            border_style="cyan",
        )
    )

    try:
        result = compare_genomes(genome1, genome2, k=kmer_size, sketch_size=sketch_size)
    except Exception as e:
        console.print(f"[bold red]☠️ Comparison failed: {e}[/bold red]")
        sys.exit(1)

    # Display results
    g1 = result["genome1"]
    g2 = result["genome2"]
    comp = result["comparison"]

    info_table = Table(title="🪨 Genome Information")
    info_table.add_column("", style="dim")
    info_table.add_column("Genome 1", style="cyan")
    info_table.add_column("Genome 2", style="green")
    info_table.add_row("Name", g1["name"], g2["name"])
    info_table.add_row("Length", f"{g1['length_bp']:,} bp", f"{g2['length_bp']:,} bp")
    info_table.add_row("Contigs", str(g1["contigs"]), str(g2["contigs"]))
    info_table.add_row("GC%", f"{g1['gc_content']}%", f"{g2['gc_content']}%")
    info_table.add_row("K-mers", f"{g1['kmers']:,}", f"{g2['kmers']:,}")
    console.print(info_table)

    result_table = Table(title="🔬 ANI Results")
    result_table.add_column("Metric", style="cyan")
    result_table.add_column("Value", style="white")
    result_table.add_row("Jaccard similarity", f"{comp['jaccard_similarity']:.6f}")
    result_table.add_row("Approximate ANI", f"{comp['approximate_ani']:.2f}%")
    result_table.add_row("Shared k-mers", f"{comp['shared_kmers']:,}")

    same = comp["same_species"]
    verdict = "[green]✅ Same species[/green]" if same else "[red]❌ Different species[/red]"
    result_table.add_row("Species call", verdict)
    result_table.add_row("Threshold", comp["threshold"])
    console.print(result_table)


@cli.command()
@click.argument("fasta", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    default=None,
    help="Output directory [default: <input_name>_domains]",
)
@click.option("-e", "--evalue", default=1e-5, type=float, help="E-value threshold (default: 1e-5)")
@click.option("--cpus", default=1, type=int, help="Number of CPUs")
@click.option(
    "--hmm-db",
    multiple=True,
    type=click.Path(exists=True, path_type=Path),
    help="HMM database(s) to search [default: auto-discover Pfam]",
)
@click.option("--width", default=70, type=int, help="Width of domain architecture map (default: 70)")
def domains(
    fasta: Path,
    output: Path | None,
    evalue: float,
    cpus: int,
    hmm_db: tuple,
    width: int,
) -> None:
    """Characterize protein domains from a FASTA file.

    Takes a protein FASTA and searches all sequences against Pfam
    HMM databases, producing a TSV table of domain hits and an
    ASCII domain architecture map.

    \b
    Example:
        darwin domains proteins.faa
        darwin domains proteins.faa -o my_results/ --cpus 4
        darwin domains proteins.faa --evalue 1e-10
    """
    from darwin.soil.nutrients import NutrientStore
    from darwin.sunlight.domains import (
        format_all_maps,
        format_tsv,
        read_protein_fasta,
        search_domains,
    )

    # Output directory
    if output is None:
        output = fasta.parent / f"{fasta.stem}_domains"
    output.mkdir(parents=True, exist_ok=True)

    # Read proteins
    console.print(
        Panel.fit(
            f"[bold green]Domain Characterization[/bold green]\n"
            f"Input: {fasta}\n"
            f"Output: {output}\n"
            f"E-value: {evalue}",
            title="🧬 Darwin Domains",
            border_style="green",
        )
    )

    proteins = read_protein_fasta(fasta)
    if not proteins:
        console.print("[red]No protein sequences found in input file.[/red]")
        sys.exit(1)

    console.print(f"  Read [bold]{len(proteins)}[/bold] protein sequences")

    # Find HMM databases
    if hmm_db:
        hmm_paths = list(hmm_db)
    else:
        soil = NutrientStore()
        available_dbs = soil.get_hmm_databases()
        hmm_paths = [db.path for db in available_dbs]

    if not hmm_paths:
        console.print(
            "[red]No HMM databases found.[/red]\n"
            "  Run [bold]darwin setup --pfam[/bold] to download Pfam-A,\n"
            "  or pass --hmm-db /path/to/database.hmm"
        )
        sys.exit(1)

    console.print(f"  Using [bold]{len(hmm_paths)}[/bold] HMM database(s): "
                  + ", ".join(p.stem for p in hmm_paths))

    # Search
    console.print("\n[bold]Searching domains...[/bold]")
    results = search_domains(proteins, hmm_paths, evalue=evalue, cpus=cpus)

    # Stats
    total_hits = sum(p.num_domains for p in results)
    proteins_with_domains = sum(1 for p in results if p.num_domains > 0)

    stats_table = Table(title="🧬 Domain Search Results")
    stats_table.add_column("Metric", style="cyan")
    stats_table.add_column("Value", style="white")
    stats_table.add_row("Total proteins", str(len(results)))
    stats_table.add_row("Proteins with domains", f"{proteins_with_domains} ({proteins_with_domains/len(results):.1%})")
    stats_table.add_row("Total domain hits", str(total_hits))
    if proteins_with_domains > 0:
        avg_domains = total_hits / proteins_with_domains
        stats_table.add_row("Avg domains per hit protein", f"{avg_domains:.1f}")
    console.print(stats_table)

    # Write TSV
    tsv_path = output / f"{fasta.stem}_domains.tsv"
    tsv_content = format_tsv(results)
    tsv_path.write_text(tsv_content)
    console.print(f"\n  📄 TSV: {tsv_path}")

    # Write domain map
    map_path = output / f"{fasta.stem}_domain_map.txt"
    map_content = format_all_maps(results, width=width)
    map_path.write_text(map_content)
    console.print(f"  🗺️  Map: {map_path}")

    # Show a preview of the first few maps
    preview_count = min(5, len(results))
    console.print(f"\n[dim]Preview (first {preview_count} proteins):[/dim]\n")
    from darwin.sunlight.domains import format_domain_map

    for prot in results[:preview_count]:
        console.print(format_domain_map(prot, width=width))

    console.print(f"\n[green bold]Done![/green bold] Results in: {output}")
