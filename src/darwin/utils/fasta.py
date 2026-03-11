"""FASTA I/O helpers using Biopython."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from darwin.models import Contig, Genome


def parse_fasta(path: Path, min_len: int = 200) -> Genome:
    """Parse a FASTA file into a Genome object, filtering short contigs."""
    contigs: list[Contig] = []
    for record in SeqIO.parse(str(path), "fasta"):
        seq = str(record.seq).upper()
        if len(seq) >= min_len:
            contigs.append(
                Contig(
                    id=record.id,
                    sequence=seq,
                    description=record.description,
                )
            )
    name = path.stem
    return Genome(name=name, contigs=contigs)


def write_fasta(contigs: list[Contig], output: Path) -> Path:
    """Write contigs to a FASTA file."""
    records = [
        SeqRecord(Seq(c.sequence), id=c.id, description=c.description)
        for c in contigs
    ]
    SeqIO.write(records, str(output), "fasta")
    return output


def write_proteins(
    features_with_seqs: list[tuple[str, str, str]],
    output: Path,
) -> Path:
    """Write protein sequences to a FASTA file.

    Args:
        features_with_seqs: list of (locus_tag, product, amino_acid_seq) tuples.
        output: path to write.
    """
    records = [
        SeqRecord(Seq(seq), id=tag, description=product)
        for tag, product, seq in features_with_seqs
    ]
    SeqIO.write(records, str(output), "fasta")
    return output
