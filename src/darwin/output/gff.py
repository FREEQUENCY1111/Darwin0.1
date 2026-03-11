"""GFF3 output format writer."""

from __future__ import annotations

from pathlib import Path

from darwin.models import Genome


def write_gff3(genome: Genome, output: Path) -> Path:
    """Write genome annotations in GFF3 format."""
    with open(output, "w") as fh:
        fh.write("##gff-version 3\n")

        for contig in genome.contigs:
            # Sequence region pragma
            fh.write(f"##sequence-region {contig.id} 1 {contig.length}\n")

            for feature in contig.features:
                # GFF3 columns: seqid source type start end score strand phase attributes
                score = f"{feature.score:.1f}" if feature.score is not None else "."
                phase = str(feature.phase) if feature.phase is not None else "."
                strand = feature.strand.value

                # Build attributes string
                attrs_parts = []
                if feature.locus_tag:
                    attrs_parts.append(f"ID={feature.locus_tag}")
                for key, val in feature.attributes.items():
                    if key == "locus_tag":
                        continue  # already in ID
                    # GFF3 escaping
                    val = val.replace(";", "%3B").replace("=", "%3D").replace(",", "%2C")
                    attrs_parts.append(f"{key}={val}")
                attrs_str = ";".join(attrs_parts)

                fh.write(
                    f"{feature.seq_id}\tDarwin\t{feature.feature_type.value}\t"
                    f"{feature.start}\t{feature.end}\t{score}\t{strand}\t{phase}\t"
                    f"{attrs_str}\n"
                )

        fh.write("###\n")

    return output
