"""
FASTA I/O — Reading and writing the raw rock material.

Parses FASTA files into Genome objects (rocks) and writes
sequences back out. This is how sunlight enters the jar.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from darwin.rocks.models import Contig, FeatureType, Genome

logger = logging.getLogger("darwin.rocks")


def parse_fasta(filepath: Path, min_length: int = 200) -> Genome:
    """
    Parse a FASTA file into a Genome.

    Supports .gz compressed files. Filters contigs below min_length.
    This is the moment sunlight enters the jar — raw sequence
    becomes structured rock.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"No rock found at: {filepath}")

    contigs: list[Contig] = []
    current_id = ""
    current_desc = ""
    current_seq: list[str] = []

    if filepath.suffix == ".gz":
        fh_ctx = gzip.open(filepath, "rt")
    else:
        fh_ctx = open(filepath)  # noqa: SIM115

    with fh_ctx as fh:
        for raw_line in fh:
            line_str: str = str(raw_line).strip()
            if not line_str:
                continue
            if line_str.startswith(">"):
                # Save previous contig
                if current_id and len("".join(current_seq)) >= min_length:
                    contigs.append(
                        Contig(
                            id=current_id,
                            sequence="".join(current_seq).upper(),
                            description=current_desc,
                        )
                    )
                # Parse header
                parts = line_str[1:].split(None, 1)
                current_id = str(parts[0])
                current_desc = str(parts[1]) if len(parts) > 1 else ""
                current_seq = []
            else:
                current_seq.append(line_str)

    # Don't forget the last contig
    if current_id and len("".join(current_seq)) >= min_length:
        contigs.append(
            Contig(
                id=current_id,
                sequence="".join(current_seq).upper(),
                description=current_desc,
            )
        )

    name = filepath.stem
    if name.endswith(".fasta") or name.endswith(".fa") or name.endswith(".fna"):
        name = Path(name).stem

    genome = Genome(name=name, contigs=contigs, source_file=filepath)
    logger.info(
        f"🪨 Parsed {genome.num_contigs} contigs, "
        f"{genome.total_length:,} bp, {genome.gc_content}% GC"
    )
    return genome


def write_fasta(genome: Genome, output: Path, wrap: int = 70) -> Path:
    """Write genome contigs to a FASTA file."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    with open(output, "w") as fh:
        for contig in genome.contigs:
            desc = f" {contig.description}" if contig.description else ""
            fh.write(f">{contig.id}{desc}\n")
            seq = contig.sequence
            for i in range(0, len(seq), wrap):
                fh.write(seq[i : i + wrap] + "\n")

    logger.info(f"📝 Wrote {genome.num_contigs} contigs to {output}")
    return output


def write_proteins(genome: Genome, output: Path, wrap: int = 70) -> Path:
    """Write translated CDS features to a protein FASTA file."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with open(output, "w") as fh:
        for feature in genome.features_of_type(FeatureType.CDS):
            if feature.translation:
                header = f">{feature.locus_tag} {feature.product}"
                fh.write(header + "\n")
                seq = feature.translation
                for i in range(0, len(seq), wrap):
                    fh.write(seq[i : i + wrap] + "\n")
                count += 1

    logger.info(f"🧬 Wrote {count} proteins to {output}")
    return output
