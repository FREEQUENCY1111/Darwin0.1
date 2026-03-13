"""
ISEScanPlant — Insertion sequence / transposon detection organism.

Feeds on: genome.loaded (raw sequence data)
Produces: mobile_elements.found (IS element features)

Like an archaeologist mapping ancient settlements — this organism
scans the genome for insertion sequences (IS elements) and
annotates them with their family, inverted repeats, and target
site duplications.

Uses ISEScan for IS element detection.
Install: conda install -c bioconda isescan
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
from darwin.water import Stream
from darwin.water.stream import Nutrient, NutrientType

logger = logging.getLogger("darwin.flora.isescan")


class ISEScanPlant(Organism):
    """ISEScan IS element detector — the archaeologist."""

    name = "isescan"
    feeds_on_nutrients = [NutrientType.GENOME_LOADED]
    produces_nutrients = [NutrientType.MOBILE_ELEMENTS_FOUND]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return self.soil.has_isescan

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """
        Run ISEScan on the genome to detect IS elements.

        Scans all contigs for insertion sequences, creates
        MOBILE_ELEMENT features, and attaches them to contigs.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        cpus = config.get("cpus", 1)
        locus_prefix = config.get("locus_tag_prefix", "DARWIN")

        if not self.can_grow():
            self.logger.warning("🏜️ isescan.py not in soil — cannot detect IS elements")
            return None

        self.logger.info(f"🏺 Scanning for IS elements in {genome.name}...")

        with tempfile.TemporaryDirectory(prefix="darwin_isescan_") as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / "input.fasta"
            output_dir = tmp / "isescan_out"

            # Write genome to temp FASTA
            with open(input_fasta, "w") as fh:
                for contig in genome.contigs:
                    fh.write(f">{contig.id}\n{contig.sequence}\n")

            # Run ISEScan
            cmd = [
                "isescan.py",
                "--seqfile", str(input_fasta),
                "--output", str(output_dir),
                "--nthread", str(cpus),
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
                    f"⚠️ ISEScan failed (exit {proc.returncode}): {err_msg}"
                )
                return Nutrient(
                    type=NutrientType.MOBILE_ELEMENTS_FOUND,
                    data={
                        "genome": genome,
                        "is_count": 0,
                        "detection_failed": True,
                        "config": config,
                    },
                    source=self.name,
                    correlation_id=nutrient.correlation_id,
                )

            # Parse ISEScan output
            is_count, families = self._parse_isescan_output(
                genome, output_dir, locus_prefix
            )

        family_str = ", ".join(sorted(families)[:5])
        if len(families) > 5:
            family_str += f" (+{len(families) - 5} more)"

        self.logger.info(
            f"🏺 Found {is_count} IS element(s) "
            f"from {len(families)} family/families"
            + (f": {family_str}" if families else "")
        )

        return Nutrient(
            type=NutrientType.MOBILE_ELEMENTS_FOUND,
            data={
                "genome": genome,
                "is_count": is_count,
                "is_families": sorted(families),
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    @staticmethod
    def _parse_isescan_output(
        genome: Genome,
        output_dir: Path,
        locus_prefix: str,
    ) -> tuple[int, set[str]]:
        """
        Parse ISEScan output and create MOBILE_ELEMENT features.

        ISEScan produces a TSV summary with columns:
          seqID, family, cluster, isBegin, isEnd, isLen,
          ncopy4is, strand, type, score, irId, irLen, tir

        Returns: (is_count, set_of_families)
        """
        contig_map = {c.id: c for c in genome.contigs}
        is_count = 0
        families: set[str] = set()

        # ISEScan output files — look for the summary TSV
        # The file is typically named <input_base>.tsv or in a subdirectory
        tsv_candidates = list(output_dir.rglob("*.tsv"))
        if not tsv_candidates:
            # Also check for .csv files
            tsv_candidates = list(output_dir.rglob("*.csv"))

        for tsv_path in tsv_candidates:
            # Skip files that are clearly not the IS summary
            if "proteome" in tsv_path.name.lower():
                continue

            try:
                with open(tsv_path, newline="") as fh:
                    # Detect delimiter
                    sample = fh.read(4096)
                    fh.seek(0)
                    dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
                    reader = csv.DictReader(fh, delimiter=dialect.delimiter)

                    for row in reader:
                        seq_id = row.get("seqID", "").strip()
                        if not seq_id or seq_id not in contig_map:
                            continue

                        family = row.get("family", "unknown").strip()
                        families.add(family)

                        try:
                            start = int(row.get("isBegin", "0"))
                            end = int(row.get("isEnd", "0"))
                        except ValueError:
                            continue

                        if start == 0 or end == 0:
                            continue

                        strand_str = row.get("strand", ".").strip()
                        strand = Strand.FORWARD if strand_str == "+" else (
                            Strand.REVERSE if strand_str == "-" else Strand.UNKNOWN
                        )

                        score_str = row.get("score", "0").strip()
                        try:
                            score = float(score_str)
                        except ValueError:
                            score = 0.0

                        is_count += 1
                        locus_tag = f"{locus_prefix}_is{is_count:04d}"

                        # Build informative note
                        notes = [f"IS family: {family}"]
                        cluster = row.get("cluster", "").strip()
                        if cluster:
                            notes.append(f"cluster: {cluster}")
                        is_type = row.get("type", "").strip()
                        if is_type:
                            notes.append(f"type: {is_type}")
                        ir_len = row.get("irLen", "").strip()
                        ir_id = row.get("irId", "").strip()
                        if ir_len and ir_len != "0":
                            notes.append(f"IR: {ir_len}bp ({ir_id}% identity)")
                        tir = row.get("tir", "").strip()
                        if tir:
                            notes.append(f"TIR: {tir}")

                        feature = Feature(
                            type=FeatureType.MOBILE_ELEMENT,
                            start=start,
                            end=end,
                            strand=strand,
                            score=score,
                            contig_id=seq_id,
                            locus_tag=locus_tag,
                            product=f"{family} family transposase",
                            inference="ab initio prediction:ISEScan",
                            note="; ".join(notes),
                        )

                        contig_map[seq_id].features.append(feature)

            except (csv.Error, KeyError) as e:
                logger.debug(f"Skipping {tsv_path.name}: {e}")
                continue

        return is_count, families
