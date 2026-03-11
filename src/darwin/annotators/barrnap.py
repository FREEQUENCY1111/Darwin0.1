"""Barrnap — ribosomal RNA prediction.

Barrnap uses Hidden Markov Models to predict 5S, 16S, and 23S rRNA
genes in prokaryotic genomes. It's fast and reliable.
"""

from __future__ import annotations

from pathlib import Path

from darwin.annotators.base import BaseAnnotator
from darwin.models import AnnotationConfig, Feature, FeatureType, Genome, Strand
from darwin.utils.runners import run_external, which_tool


class BarrnapAnnotator(BaseAnnotator):
    name = "barrnap"

    def check_dependencies(self) -> bool:
        try:
            which_tool("barrnap")
            return True
        except Exception:
            return False

    def run(self, genome: Genome) -> Genome:
        self.log.info("[bold blue]Running Barrnap rRNA prediction...[/]")

        input_fasta = self._write_input(genome)
        gff_out = self.work_dir / "barrnap.gff3"

        kingdom = self.config.kingdom
        if kingdom not in ("bac", "arc", "euk", "mito"):
            kingdom = "bac"

        cmd = [
            "barrnap",
            "--kingdom", kingdom,
            "--threads", str(self.config.cpus),
            "--reject", "0.2",  # minimum length proportion
            "--quiet",
            str(input_fasta),
        ]
        result = run_external(cmd, capture=True)

        # Barrnap writes GFF3 to stdout
        gff_out.write_text(result.stdout)
        genome = self._parse_gff(genome, gff_out)

        rrna_count = sum(
            1
            for c in genome.contigs
            for f in c.features
            if f.feature_type == FeatureType.RRNA
        )
        self.log.info(f"  Found [green]{rrna_count}[/] rRNA genes")
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

    def _parse_gff(self, genome: Genome, gff_path: Path) -> Genome:
        """Parse Barrnap GFF3 output."""
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
                start = int(parts[3])
                end = int(parts[4])
                score = float(parts[5]) if parts[5] != "." else None
                strand = Strand.FORWARD if parts[6] == "+" else Strand.REVERSE

                # Parse attributes for product name
                attrs_raw = {}
                for item in parts[8].split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        attrs_raw[k] = v

                # Barrnap names like "16S_rRNA", "23S_rRNA", "5S_rRNA"
                product = attrs_raw.get("Name", attrs_raw.get("product", "rRNA"))
                # Clean up: "16S_rRNA" -> "16S ribosomal RNA"
                product_clean = product.replace("_rRNA", " ribosomal RNA").replace("_", " ")

                gene_counter += 1
                locus_tag = f"{self.config.locus_tag_prefix}_r{gene_counter:04d}"

                feature = Feature(
                    seq_id=seq_id,
                    feature_type=FeatureType.RRNA,
                    start=start,
                    end=end,
                    strand=strand,
                    score=score,
                    attributes={
                        "locus_tag": locus_tag,
                        "product": product_clean,
                        "inference": "profile:Barrnap",
                    },
                )

                if seq_id in contig_map:
                    contig_map[seq_id].features.append(feature)

        return genome
