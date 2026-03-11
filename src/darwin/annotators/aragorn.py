"""Aragorn — tRNA and tmRNA detection.

Aragorn detects tRNA genes and tmRNA (transfer-messenger RNA) in
nucleotide sequences. It's faster and often more accurate than
tRNAscan-SE for prokaryotic genomes.
"""

from __future__ import annotations

import re
from pathlib import Path

from darwin.annotators.base import BaseAnnotator
from darwin.models import AnnotationConfig, Feature, FeatureType, Genome, Strand
from darwin.utils.runners import run_external, which_tool


class AragornAnnotator(BaseAnnotator):
    name = "aragorn"

    def check_dependencies(self) -> bool:
        try:
            which_tool("aragorn")
            return True
        except Exception:
            return False

    def run(self, genome: Genome) -> Genome:
        self.log.info("[bold blue]Running Aragorn tRNA/tmRNA detection...[/]")

        input_fasta = self._write_input(genome)
        output_file = self.work_dir / "aragorn.txt"

        cmd = [
            "aragorn",
            "-t",                    # search for tRNA
            "-m",                    # search for tmRNA
            "-gc" + str(self.config.translation_table),
            "-l",                    # assume linear topology
            "-w",                    # batch mode output
            "-o", str(output_file),
            str(input_fasta),
        ]
        run_external(cmd)

        genome = self._parse_output(genome, output_file)

        trna_count = sum(
            1
            for c in genome.contigs
            for f in c.features
            if f.feature_type == FeatureType.TRNA
        )
        tmrna_count = sum(
            1
            for c in genome.contigs
            for f in c.features
            if f.feature_type == FeatureType.TMRNA
        )
        self.log.info(f"  Found [green]{trna_count}[/] tRNAs, [green]{tmrna_count}[/] tmRNAs")
        return genome

    def _write_input(self, genome: Genome) -> Path:
        fasta_path = self.work_dir / "input.fasta"
        with open(fasta_path, "w") as fh:
            for contig in genome.contigs:
                fh.write(f">{contig.id}\n")
                seq = contig.sequence
                for i in range(0, len(seq), 80):
                    fh.write(seq[i : i + 80] + "\n")
        return fasta_path

    def _parse_output(self, genome: Genome, output_path: Path) -> Genome:
        """Parse Aragorn's text output format."""
        contig_map = {c.id: c for c in genome.contigs}
        current_contig: str | None = None
        gene_counter = 0

        # Aragorn output format:
        # >contig_id
        # N genes found
        # N    tRNA-Xxx          c[start,end]     ...

        with open(output_path) as fh:
            for line in fh:
                line = line.strip()

                # New contig header
                if line.startswith(">"):
                    current_contig = line[1:].split()[0]
                    continue

                if current_contig is None:
                    continue

                # Match tRNA/tmRNA lines
                # Examples:
                #   1    tRNA-Ala           c[1234,1305]      34      (tgc)
                #   2    tmRNA              [5678,5901]       ...
                match = re.match(
                    r"\s*\d+\s+(tRNA-\w+|tmRNA)\s+(c?)\[(\d+),(\d+)\]",
                    line,
                )
                if not match:
                    continue

                rna_type = match.group(1)
                complement = match.group(2)
                start = int(match.group(3))
                end = int(match.group(4))
                strand = Strand.REVERSE if complement == "c" else Strand.FORWARD

                # Determine feature type and product
                if rna_type.startswith("tRNA"):
                    feat_type = FeatureType.TRNA
                    amino_acid = rna_type.replace("tRNA-", "")
                    product = f"tRNA-{amino_acid}"

                    # Extract anticodon if present
                    anticodon_match = re.search(r"\((\w{3})\)", line)
                    anticodon = anticodon_match.group(1) if anticodon_match else ""
                else:
                    feat_type = FeatureType.TMRNA
                    product = "tmRNA"
                    anticodon = ""

                gene_counter += 1
                locus_tag = f"{self.config.locus_tag_prefix}_t{gene_counter:04d}"

                attrs = {
                    "locus_tag": locus_tag,
                    "product": product,
                    "inference": "profile:Aragorn",
                }
                if anticodon:
                    attrs["anticodon"] = anticodon

                feature = Feature(
                    seq_id=current_contig,
                    feature_type=feat_type,
                    start=start,
                    end=end,
                    strand=strand,
                    attributes=attrs,
                )

                if current_contig in contig_map:
                    contig_map[current_contig].features.append(feature)

        return genome
