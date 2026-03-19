"""TSV feature table output format writer."""

from __future__ import annotations

import csv
from pathlib import Path

from darwin.rocks.models import Genome


def write_tsv(genome: Genome, output: Path) -> Path:
    """Write genome annotations as a tab-separated feature table."""
    headers = [
        "locus_tag", "type", "contig", "start", "end", "strand",
        "length", "product", "score", "inference", "db_xref",
        "gene", "operon", "go_terms", "ipr_ids", "structure_hit", "note",
    ]

    with open(output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(headers)

        for contig in genome.contigs:
            sorted_features = sorted(contig.features, key=lambda f: f.start)
            for f in sorted_features:
                # Extract operon info from note
                operon = ""
                note_parts = []
                if f.note:
                    for part in f.note.split(";"):
                        if part.startswith("operon="):
                            operon = part.split("=", 1)[1]
                        elif part.startswith("operon_pos="):
                            operon += f":{part.split('=', 1)[1]}"
                        else:
                            note_parts.append(part)

                writer.writerow([
                    f.locus_tag,
                    f.type.value,
                    contig.id,
                    f.start,
                    f.end,
                    f.strand.value,
                    f.length,
                    f.product,
                    f"{f.score:.1f}" if f.score else "",
                    f.inference,
                    ";".join(f.db_xref) if f.db_xref else "",
                    f.gene,
                    operon,
                    ";".join(f.go_terms) if f.go_terms else "",
                    ";".join(f.ipr_ids) if f.ipr_ids else "",
                    f.structure_hit,
                    ";".join(note_parts),
                ])

    return output
