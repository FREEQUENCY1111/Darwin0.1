"""
ABRicatePlant — Antimicrobial resistance gene detection.

Feeds on: genome.loaded (raw sequence data)
Produces: resistance_genes.found (AMR gene annotations)

Like antibiotic-resistant bacteria in a jar — they carry resistance
genes that reveal what chemical warfare they've survived. ABRicate
screens the genome against curated resistance databases to find
these molecular shields.

Uses ABRicate for multi-database resistance gene screening.
Install: conda install -c bioconda abricate
"""

from __future__ import annotations

import asyncio
import csv
import logging
import tempfile
from pathlib import Path

from darwin.flora.base import Organism
from darwin.rocks.models import Feature, FeatureType, Genome, Strand
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.flora.abricate")


class ABRicatePlant(Organism):
    """ABRicate AMR screening — molecular resistance detection."""

    name = "abricate"
    feeds_on_nutrients = [NutrientType.GENOME_LOADED]
    produces_nutrients = [NutrientType.RESISTANCE_GENES_FOUND]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        """Check if ABRicate is available in soil."""
        return self.soil.has_abricate

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """
        Screen genome for antimicrobial resistance genes.

        Runs ABRicate with the NCBI AMRFinder database by default.
        Creates AMR_GENE features for each hit with resistance class,
        identity, and database cross-references.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        cpus = config.get("cpus", 1)
        locus_prefix = config.get("locus_tag_prefix", "DARWIN")

        if not self.can_grow():
            self.logger.warning("🏜️ ABRicate not in soil — cannot detect AMR genes")
            return Nutrient(
                type=NutrientType.RESISTANCE_GENES_FOUND,
                data={
                    "genome": genome,
                    "amr_count": 0,
                    "resistance_classes": [],
                    "detection_skipped": True,
                    "config": config,
                },
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        self.logger.info(f"💊 Screening {genome.name} for antimicrobial resistance genes...")

        with tempfile.TemporaryDirectory(prefix="darwin_abricate_") as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / "input.fasta"
            output_tsv = tmp / "abricate_results.tsv"

            # Write genome to temp FASTA
            with open(input_fasta, "w") as fh:
                for contig in genome.contigs:
                    fh.write(f">{contig.id}\n{contig.sequence}\n")

            # Run ABRicate
            cmd = [
                "abricate",
                "--quiet",
                "--threads", str(cpus),
                "--minid", "80",
                "--mincov", "60",
                str(input_fasta),
            ]

            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, stderr = await proc.communicate()

            if proc.returncode != 0:
                err_msg = stderr.decode()[:500]
                self.logger.warning(
                    f"⚠️ ABRicate failed (exit {proc.returncode}): {err_msg}"
                )
                return Nutrient(
                    type=NutrientType.RESISTANCE_GENES_FOUND,
                    data={
                        "genome": genome,
                        "amr_count": 0,
                        "resistance_classes": [],
                        "detection_failed": True,
                        "config": config,
                    },
                    source=self.name,
                    correlation_id=nutrient.correlation_id,
                )

            # ABRicate writes to stdout as TSV
            with open(output_tsv, "wb") as fh:
                fh.write(stdout)

            amr_count, resistance_classes = self._parse_abricate_output(
                genome, output_tsv, locus_prefix
            )

        self.logger.info(
            f"💊 Found {amr_count} AMR gene(s) across "
            f"{len(resistance_classes)} resistance class(es)"
        )

        return Nutrient(
            type=NutrientType.RESISTANCE_GENES_FOUND,
            data={
                "genome": genome,
                "amr_count": amr_count,
                "resistance_classes": sorted(resistance_classes),
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    @staticmethod
    def _parse_abricate_output(
        genome: Genome,
        output_path: Path,
        locus_prefix: str,
    ) -> tuple[int, set[str]]:
        """
        Parse ABRicate TSV output and create AMR_GENE features.

        ABRicate output columns (tab-delimited):
        #FILE  SEQUENCE  START  END  STRAND  GENE  COVERAGE  COVERAGE_MAP
        GAPS  %COVERAGE  %IDENTITY  DATABASE  ACCESSION  PRODUCT  RESISTANCE

        Returns (amr_count, resistance_classes_set).
        """
        contig_map = {c.id: c for c in genome.contigs}
        amr_count = 0
        resistance_classes: set[str] = set()

        if not output_path.exists():
            return 0, resistance_classes

        try:
            with open(output_path, newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    seq_id = row.get("SEQUENCE", "").strip()
                    if not seq_id or seq_id not in contig_map:
                        continue

                    try:
                        start = int(row.get("START", "0"))
                        end = int(row.get("END", "0"))
                    except ValueError:
                        continue

                    if start == 0 or end == 0:
                        continue

                    amr_count += 1
                    locus_tag = f"{locus_prefix}_amr{amr_count:04d}"

                    gene_name = row.get("GENE", "unknown").strip()
                    product = row.get("PRODUCT", gene_name).strip()
                    resistance = row.get("RESISTANCE", "").strip()
                    identity = row.get("%IDENTITY", "").strip()
                    coverage = row.get("%COVERAGE", "").strip()
                    database = row.get("DATABASE", "").strip()
                    accession = row.get("ACCESSION", "").strip()

                    # Track resistance classes
                    if resistance and resistance != "-":
                        for cls in resistance.split(";"):
                            resistance_classes.add(cls.strip())

                    # Determine strand
                    strand_str = row.get("STRAND", "+").strip()
                    strand = Strand.REVERSE if strand_str == "-" else Strand.FORWARD

                    # Build note with metadata
                    note_parts = []
                    if resistance and resistance != "-":
                        note_parts.append(f"resistance class: {resistance}")
                    if identity:
                        note_parts.append(f"identity: {identity}%")
                    if coverage:
                        note_parts.append(f"coverage: {coverage}%")
                    if database:
                        note_parts.append(f"database: {database}")

                    # Build db_xref
                    db_xref = []
                    if accession and accession != "-":
                        db_xref.append(f"AMR:{accession}")

                    feature = Feature(
                        type=FeatureType.AMR_GENE,
                        start=start,
                        end=end,
                        strand=strand,
                        contig_id=seq_id,
                        locus_tag=locus_tag,
                        product=product or gene_name,
                        gene=gene_name,
                        inference=f"similar to AA sequence:ABRicate:{accession}",
                        note="; ".join(note_parts),
                        db_xref=db_xref,
                    )
                    contig_map[seq_id].features.append(feature)

        except (csv.Error, KeyError, ValueError) as e:
            logger.debug(f"Error parsing ABRicate output: {e}")

        return amr_count, resistance_classes
