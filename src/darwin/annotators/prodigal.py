"""Prodigal — ab initio prokaryotic gene prediction.

Prodigal predicts protein-coding genes (CDS) using dynamic programming
on the input contigs. It's the standard first step in any prokaryotic
annotation pipeline.

Modes:
  - single: for complete/draft genomes (default)
  - meta:   for metagenomic fragments
"""

from __future__ import annotations

import re
from pathlib import Path

from darwin.annotators.base import BaseAnnotator
from darwin.models import AnnotationConfig, Feature, FeatureType, Genome, Strand
from darwin.utils.runners import run_external, which_tool


class ProdigalAnnotator(BaseAnnotator):
    name = "prodigal"

    def check_dependencies(self) -> bool:
        try:
            which_tool("prodigal")
            return True
        except Exception:
            return False

    def run(self, genome: Genome) -> Genome:
        self.log.info("[bold blue]Running Prodigal gene prediction...[/]")

        # Write genome contigs to temp FASTA
        input_fasta = self._write_input(genome)
        gff_out = self.work_dir / "prodigal.gff"
        proteins_out = self.work_dir / "prodigal.faa"
        genes_out = self.work_dir / "prodigal.fna"

        # Build command
        mode = "meta" if self.config.metagenome else "single"
        cmd = [
            "prodigal",
            "-i", str(input_fasta),
            "-o", str(gff_out),
            "-f", "gff",          # GFF output
            "-a", str(proteins_out),  # protein translations
            "-d", str(genes_out),     # nucleotide gene sequences
            "-p", mode,
            "-g", str(self.config.translation_table),
            "-q",  # quiet
        ]
        run_external(cmd)

        # Parse GFF and add features
        genome = self._parse_gff(genome, gff_out)

        cds_count = sum(
            1
            for c in genome.contigs
            for f in c.features
            if f.feature_type == FeatureType.CDS
        )
        self.log.info(f"  Found [green]{cds_count}[/] protein-coding genes")
        return genome

    def _write_input(self, genome: Genome) -> Path:
        """Write genome contigs to a FASTA for Prodigal."""
        fasta_path = self.work_dir / "input.fasta"
        with open(fasta_path, "w") as fh:
            for contig in genome.contigs:
                fh.write(f">{contig.id} {contig.description}\n")
                # Write sequence in 80-char lines
                seq = contig.sequence
                for i in range(0, len(seq), 80):
                    fh.write(seq[i : i + 80] + "\n")
        return fasta_path

    def _parse_gff(self, genome: Genome, gff_path: Path) -> Genome:
        """Parse Prodigal GFF output and attach CDS features to contigs."""
        # Build a lookup for contigs
        contig_map = {c.id: c for c in genome.contigs}
        gene_counter = 0

        with open(gff_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                seq_id = parts[0]
                feat_type = parts[2]
                if feat_type != "CDS":
                    continue

                start = int(parts[3])
                end = int(parts[4])
                score = float(parts[5]) if parts[5] != "." else None
                strand = Strand.FORWARD if parts[6] == "+" else Strand.REVERSE
                phase = int(parts[7]) if parts[7] != "." else 0

                # Parse attributes
                attrs = self._parse_gff_attributes(parts[8])
                gene_counter += 1
                locus_tag = f"{self.config.locus_tag_prefix}_{gene_counter:05d}"
                attrs["locus_tag"] = locus_tag
                attrs["inference"] = "ab initio prediction:Prodigal"

                feature = Feature(
                    seq_id=seq_id,
                    feature_type=FeatureType.CDS,
                    start=start,
                    end=end,
                    strand=strand,
                    score=score,
                    phase=phase,
                    attributes=attrs,
                )

                if seq_id in contig_map:
                    contig_map[seq_id].features.append(feature)

        return genome

    @staticmethod
    def _parse_gff_attributes(attr_str: str) -> dict[str, str]:
        """Parse GFF attribute column (key=value;key=value)."""
        attrs: dict[str, str] = {}
        for item in attr_str.split(";"):
            item = item.strip()
            if "=" in item:
                key, val = item.split("=", 1)
                attrs[key] = val
        return attrs
