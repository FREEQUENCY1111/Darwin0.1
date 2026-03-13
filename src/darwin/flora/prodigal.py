"""
Prodigal Plant — Gene calling organism.

Feeds on: genome.loaded (raw sequence data)
Produces: genes.called (predicted coding sequences)

Like the first plant in a new ecosystem — it's the primary
producer that converts raw sunlight into organic matter
that everything else feeds on.
"""

from __future__ import annotations

import asyncio
import logging
import tempfile
from pathlib import Path

from darwin.flora.base import Organism
from darwin.rocks.models import Feature, FeatureType, Genome, Strand
from darwin.soil.nutrients import NutrientStore
from darwin.water import Stream
from darwin.water.stream import Nutrient, NutrientType

logger = logging.getLogger("darwin.flora.prodigal")


class ProdigalPlant(Organism):
    """Prodigal gene caller — the primary producer."""

    name = "prodigal"
    feeds_on_nutrients = [NutrientType.GENOME_LOADED]
    produces_nutrients = [NutrientType.GENES_CALLED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return self.soil.has_prodigal

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """
        Run Prodigal on the genome.

        Takes raw contigs, predicts genes, attaches CDS features
        to the genome's contigs, then releases genes.called.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        metagenome = config.get("metagenome_mode", False)
        translation_table = config.get("translation_table", 11)
        locus_prefix = config.get("locus_tag_prefix", "DARWIN")

        if not self.can_grow():
            self.logger.warning("🏜️ Prodigal not in soil — cannot grow")
            return None

        self.logger.info(f"🌱 Growing genes from {genome.name}...")

        with tempfile.TemporaryDirectory(prefix="darwin_prodigal_") as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / "input.fasta"
            output_gff = tmp / "output.gff"
            output_proteins = tmp / "proteins.faa"

            # Write genome to temp FASTA
            with open(input_fasta, "w") as fh:
                for contig in genome.contigs:
                    fh.write(f">{contig.id}\n{contig.sequence}\n")

            # Run prodigal
            cmd = [
                "prodigal",
                "-i",
                str(input_fasta),
                "-o",
                str(output_gff),
                "-f",
                "gff",
                "-a",
                str(output_proteins),
                "-g",
                str(translation_table),
            ]
            if metagenome:
                cmd.extend(["-p", "meta"])

            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, stderr = await proc.communicate()

            if proc.returncode != 0:
                raise RuntimeError(f"Prodigal failed (exit {proc.returncode}): {stderr.decode()}")

            # Parse GFF output
            features = self._parse_gff(output_gff, locus_prefix)

            # Parse protein translations
            translations = self._parse_proteins(output_proteins)

            # Attach translations to features by index
            # (Prodigal guarantees same ordering in GFF and protein FASTA)
            for i, feature in enumerate(features):
                if i < len(translations):
                    feature.translation = translations[i]

            # Attach features to contigs
            gene_count = 0
            for feature in features:
                for contig in genome.contigs:
                    if contig.id == feature.contig_id:
                        contig.features.append(feature)
                        gene_count += 1
                        break

        self.logger.info(f"🌿 Grew {gene_count} genes from {genome.name}")

        return Nutrient(
            type=NutrientType.GENES_CALLED,
            data={
                "genome": genome,
                "gene_count": gene_count,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    def _parse_gff(self, gff_path: Path, locus_prefix: str) -> list[Feature]:
        """Parse Prodigal GFF3 output into Feature objects."""
        features = []
        gene_num = 0

        with open(gff_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue

                contig_id = parts[0]
                ftype = parts[2]
                if ftype != "CDS":
                    continue

                gene_num += 1
                start = int(parts[3])
                end = int(parts[4])
                score = float(parts[5]) if parts[5] != "." else 0.0
                strand = Strand.FORWARD if parts[6] == "+" else Strand.REVERSE
                locus_tag = f"{locus_prefix}_{gene_num:05d}"

                features.append(
                    Feature(
                        type=FeatureType.CDS,
                        start=start,
                        end=end,
                        strand=strand,
                        score=score,
                        contig_id=contig_id,
                        locus_tag=locus_tag,
                        product="hypothetical protein",
                        inference="ab initio prediction:Prodigal",
                    )
                )

        return features

    def _parse_proteins(self, protein_path: Path) -> list[str]:
        """Parse Prodigal protein FASTA, return translations in GFF order."""
        if not protein_path.exists():
            return []

        current_seq: list[str] = []

        # Prodigal headers look like: >contig_1_1 # 1 # 200 # 1 # ...
        # We need to map these to our locus tags by position
        # For simplicity, we'll use a counter that matches GFF parsing order
        seqs: list[str] = []

        with open(protein_path) as fh:
            for line in fh:
                if line.startswith(">"):
                    if current_seq:
                        seqs.append("".join(current_seq).rstrip("*"))
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_seq:
                seqs.append("".join(current_seq).rstrip("*"))

        # Map to locus tags by index
        # (assumes same ordering as GFF, which Prodigal guarantees)
        return seqs
