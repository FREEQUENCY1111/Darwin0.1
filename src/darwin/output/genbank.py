"""GenBank flat file output format writer."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from darwin.models import FeatureType, Genome


def write_genbank(genome: Genome, output: Path) -> Path:
    """Write genome annotations in GenBank format using Biopython."""
    records: list[SeqRecord] = []

    for contig in genome.contigs:
        record = SeqRecord(
            seq=Seq(contig.sequence),
            id=contig.id,
            name=contig.id[:16],  # GenBank name limit
            description=f"{genome.name} {contig.id}",
            annotations={
                "molecule_type": "DNA",
                "topology": "linear",
                "data_file_division": "BCT",
                "date": datetime.now().strftime("%d-%b-%Y").upper(),
                "source": genome.name,
                "organism": genome.name,
            },
        )

        for feature in contig.features:
            # Convert strand
            strand = 1 if feature.strand.value == "+" else -1

            location = FeatureLocation(
                feature.start - 1,  # Biopython uses 0-based
                feature.end,
                strand=strand,
            )

            qualifiers: dict[str, list[str]] = {}
            if feature.locus_tag:
                qualifiers["locus_tag"] = [feature.locus_tag]
            if feature.product:
                qualifiers["product"] = [feature.product]
            if "inference" in feature.attributes:
                qualifiers["inference"] = [feature.attributes["inference"]]
            if "anticodon" in feature.attributes:
                qualifiers["anticodon"] = [feature.attributes["anticodon"]]

            # Translate CDS features
            if feature.feature_type == FeatureType.CDS:
                qualifiers["codon_start"] = [str((feature.phase or 0) + 1)]
                qualifiers["transl_table"] = ["11"]
                try:
                    nuc = contig.sequence[feature.start - 1 : feature.end]
                    if strand == -1:
                        nuc = str(Seq(nuc).reverse_complement())
                    protein = str(Seq(nuc).translate(table=11, to_stop=True))
                    qualifiers["translation"] = [protein]
                except Exception:
                    pass

            bio_feature = SeqFeature(
                location=location,
                type=feature.feature_type.value,
                qualifiers=qualifiers,
            )
            record.features.append(bio_feature)

        records.append(record)

    SeqIO.write(records, str(output), "genbank")
    return output
