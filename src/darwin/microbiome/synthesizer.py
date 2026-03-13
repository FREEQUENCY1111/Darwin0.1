"""
Synthesizer — Output generation decomposer.

Feeds on: annotation.ready
Produces: output.written → triggers EQUILIBRIUM

Like the final stage of decomposition where complex organic
matter becomes stable humus — the Synthesizer takes all the
raw annotation data and crystallizes it into final output files.

When it finishes, equilibrium is reached. The jar is complete.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

from darwin.flora.base import Organism
from darwin.output.tsv import write_tsv
from darwin.rocks.fasta import write_fasta, write_proteins
from darwin.rocks.models import FeatureType, Genome
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.microbiome.synthesizer")


class Synthesizer(Organism):
    """
    Output synthesizer — crystallizes results into files.

    Generates: GFF3, GenBank, TSV, JSON summary, protein FASTA,
    and the audit trail.
    """

    name = "synthesizer"
    feeds_on_nutrients = [NutrientType.ANNOTATION_READY]
    produces_nutrients = [NutrientType.OUTPUT_WRITTEN]

    def __init__(
        self, stream: Stream, soil: NutrientStore, output_dir: Path = Path("darwin_output")
    ) -> None:
        super().__init__(stream, soil)
        self.output_dir = Path(output_dir)

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """Generate all output files."""
        genome: Genome = nutrient.data["genome"]
        enrichment = nutrient.data.get("enrichment", {})
        qc = nutrient.data.get("qc", [])
        config = nutrient.data.get("config", {})

        output_dir = Path(config.get("output_dir", self.output_dir))
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"📋 Synthesizing output to {output_dir}...")

        files_written = []

        # 1. GFF3
        gff_path = output_dir / f"{genome.name}.gff3"
        self._write_gff3(genome, gff_path)
        files_written.append(str(gff_path))

        # 2. GenBank
        gbk_path = output_dir / f"{genome.name}.gbk"
        self._write_genbank(genome, gbk_path, config)
        files_written.append(str(gbk_path))

        # 3. TSV feature table
        tsv_path = output_dir / f"{genome.name}_features.tsv"
        write_tsv(genome, tsv_path)
        files_written.append(str(tsv_path))

        # 4. Protein FASTA
        proteins_path = output_dir / f"{genome.name}_proteins.faa"
        write_proteins(genome, proteins_path)
        files_written.append(str(proteins_path))

        # 5. Nucleotide FASTA (clean copy)
        fasta_path = output_dir / f"{genome.name}.fna"
        write_fasta(genome, fasta_path)
        files_written.append(str(fasta_path))

        # 6. JSON summary
        json_path = output_dir / f"{genome.name}_summary.json"
        summary = {
            "genome": genome.summary(),
            "enrichment": enrichment,
            "qc_checks": qc,
            "files": files_written,
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "water_cycle": self.stream.get_sediment_summary(),
        }
        with open(json_path, "w") as fh:
            json.dump(summary, fh, indent=2, default=str)
        files_written.append(str(json_path))

        self.logger.info(f"📋 Wrote {len(files_written)} output files")

        return Nutrient(
            type=NutrientType.OUTPUT_WRITTEN,
            data={
                "genome": genome,
                "files": files_written,
                "output_dir": str(output_dir),
                "summary": summary,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    def _write_gff3(self, genome: Genome, path: Path) -> None:
        """Write GFF3 format output with Dbxref and species pragma."""
        with open(path, "w") as fh:
            fh.write("##gff-version 3\n")

            # Species pragma if taxonomy inferred
            if genome.taxonomy:
                fh.write(f"##species {genome.taxonomy}\n")

            for contig in genome.contigs:
                fh.write(f"##sequence-region {contig.id} 1 {contig.length}\n")
                sorted_features = sorted(contig.features, key=lambda f: f.start)
                for f in sorted_features:
                    strand = f.strand.value
                    attrs = f"ID={f.locus_tag};product={f.product}"
                    if f.inference:
                        attrs += f";inference={f.inference}"
                    if f.gene:
                        attrs += f";gene={f.gene}"
                    if f.db_xref:
                        attrs += f";Dbxref={','.join(f.db_xref)}"
                    if f.note:
                        attrs += f";note={f.note}"
                    fh.write(
                        f"{contig.id}\tDarwin\t{f.type.value}\t"
                        f"{f.start}\t{f.end}\t{f.score:.1f}\t"
                        f"{strand}\t.\t{attrs}\n"
                    )
        logger.info(f"  Wrote GFF3: {path}")

    def _write_genbank(self, genome: Genome, path: Path, config: dict | None = None) -> None:
        """Write GenBank flat file format output with db_xref and signal_peptide."""
        transl_table = config.get("translation_table", 11) if config else 11

        with open(path, "w") as fh:
            for contig in genome.contigs:
                # Header — set topology based on replicon circularity
                topology = "circular" if contig.is_circular else "linear"
                division = "BCT"
                fh.write(
                    f"LOCUS       {contig.id:<16} {contig.length} bp    DNA     {topology:<8} {division}\n"
                )
                fh.write(f"DEFINITION  {genome.name} {contig.description}\n")
                fh.write(f"ACCESSION   {contig.id}\n")
                fh.write(f"VERSION     {contig.id}\n")
                fh.write(f"SOURCE      {genome.organism or 'Unknown organism'}\n")
                fh.write(f"  ORGANISM  {genome.organism or 'Unknown organism'}\n")
                if genome.taxonomy:
                    fh.write(f"            {genome.taxonomy}\n")
                fh.write("FEATURES             Location/Qualifiers\n")
                fh.write(f"     source          1..{contig.length}\n")
                fh.write(f'                     /organism="{genome.organism or "Unknown"}"\n')
                fh.write('                     /mol_type="genomic DNA"\n')
                if contig.replicon_type == "plasmid":
                    fh.write('                     /plasmid="unnamed"\n')
                if contig.rep_type:
                    fh.write(f'                     /note="replicon type: {contig.rep_type}"\n')

                sorted_features = sorted(contig.features, key=lambda f: f.start)
                for f in sorted_features:
                    loc = f.location_str
                    fh.write(f"     {f.type.value:<16}{loc}\n")
                    if f.locus_tag:
                        fh.write(f'                     /locus_tag="{f.locus_tag}"\n')
                    fh.write(f'                     /product="{f.product}"\n')
                    if f.inference:
                        fh.write(f'                     /inference="{f.inference}"\n')
                    if f.gene:
                        fh.write(f'                     /gene="{f.gene}"\n')
                    # db_xref qualifiers
                    for xref in f.db_xref:
                        fh.write(f'                     /db_xref="{xref}"\n')
                    if f.note:
                        fh.write(f'                     /note="{f.note}"\n')
                    if f.translation and f.type == FeatureType.CDS:
                        fh.write(f'                     /transl_table={transl_table}\n')
                        # Wrap translation at 58 chars per line
                        fh.write('                     /translation="')
                        seq = f.translation
                        fh.write(seq[:58])
                        for i in range(58, len(seq), 58):
                            fh.write(f"\n                     {seq[i : i + 58]}")
                        fh.write('"\n')

                # Sequence
                fh.write("ORIGIN\n")
                seq = contig.sequence.lower()
                for i in range(0, len(seq), 60):
                    chunk = seq[i : i + 60]
                    parts = [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
                    fh.write(f"{i + 1:>9} {' '.join(parts)}\n")
                fh.write("//\n")

        logger.info(f"  Wrote GenBank: {path}")
