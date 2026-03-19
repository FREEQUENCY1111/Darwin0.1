"""
InterProPlant — Protein functional annotation via InterProScan.

Feeds on: proteins.found (annotated proteins with translations)
Produces: functions.annotated (GO terms, IPR IDs, pathway annotations)

Like a mycorrhizal network in the jar — it connects each protein
to the vast InterPro consortium of databases (Pfam, TIGRFAM, CDD,
Gene3D, SMART, PANTHER, etc.), enriching each gene with functional
context: what it does (GO terms), what family it belongs to (IPR),
and what pathways it participates in.

Uses InterProScan for multi-database functional annotation.
Install: download from https://interproscan-docs.readthedocs.io/
"""

from __future__ import annotations

import asyncio
import csv
import logging
import tempfile
from pathlib import Path

from darwin.flora.base import Organism
from darwin.rocks.models import FeatureType, Genome
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.flora.interpro")


class InterProPlant(Organism):
    """InterProScan functional annotator — the mycorrhizal network."""

    name = "interproscan"
    feeds_on_nutrients = [NutrientType.PROTEINS_FOUND]
    produces_nutrients = [NutrientType.FUNCTIONS_ANNOTATED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        """Check if InterProScan is available in soil."""
        return self.soil.has_interproscan

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """
        Run InterProScan on all protein sequences.

        InterProScan scans proteins against multiple member databases
        and assigns GO terms, IPR identifiers, and pathway annotations.
        Output is parsed from TSV format.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        cpus = config.get("cpus", 1)

        if not self.can_grow():
            self.logger.warning(
                "🏜️ InterProScan not in soil — cannot annotate protein functions"
            )
            return Nutrient(
                type=NutrientType.FUNCTIONS_ANNOTATED,
                data={
                    "genome": genome,
                    "go_count": 0,
                    "ipr_count": 0,
                    "annotated_count": 0,
                    "detection_skipped": True,
                    "config": config,
                },
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        self.logger.info(
            f"🍄 Running InterProScan on {genome.name} proteins..."
        )

        cds_features = genome.features_of_type(FeatureType.CDS)
        proteins = {f.locus_tag: f for f in cds_features if f.translation}

        if not proteins:
            self.logger.warning("No proteins to scan")
            return Nutrient(
                type=NutrientType.FUNCTIONS_ANNOTATED,
                data={
                    "genome": genome,
                    "go_count": 0,
                    "ipr_count": 0,
                    "annotated_count": 0,
                    "config": config,
                },
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        with tempfile.TemporaryDirectory(prefix="darwin_iprscan_") as tmpdir:
            tmp = Path(tmpdir)
            input_faa = tmp / "proteins.faa"
            output_base = tmp / "iprscan_results"

            # Write protein FASTA
            with open(input_faa, "w") as fh:
                for tag, feat in proteins.items():
                    fh.write(f">{tag}\n{feat.translation}\n")

            # Run InterProScan
            cmd = [
                "interproscan.sh",
                "-i", str(input_faa),
                "-o", str(output_base) + ".tsv",
                "-f", "TSV",
                "--goterms",
                "--iprlookup",
                "--pathways",
                "--cpu", str(cpus),
                "--disable-precalc",
            ]

            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, stderr = await proc.communicate()

            if proc.returncode != 0:
                err_msg = stderr.decode()[:500] or stdout.decode()[:500]
                self.logger.warning(
                    f"⚠️ InterProScan failed (exit {proc.returncode}): {err_msg}"
                )
                return Nutrient(
                    type=NutrientType.FUNCTIONS_ANNOTATED,
                    data={
                        "genome": genome,
                        "go_count": 0,
                        "ipr_count": 0,
                        "annotated_count": 0,
                        "detection_failed": True,
                        "config": config,
                    },
                    source=self.name,
                    correlation_id=nutrient.correlation_id,
                )

            go_count, ipr_count, annotated_count = self._parse_iprscan_output(
                proteins, output_base.with_suffix(".tsv")
            )

        self.logger.info(
            f"🍄 InterProScan: {annotated_count} proteins annotated, "
            f"{go_count} GO terms, {ipr_count} IPR IDs assigned"
        )

        return Nutrient(
            type=NutrientType.FUNCTIONS_ANNOTATED,
            data={
                "genome": genome,
                "go_count": go_count,
                "ipr_count": ipr_count,
                "annotated_count": annotated_count,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    @staticmethod
    def _parse_iprscan_output(
        proteins: dict[str, object],
        tsv_path: Path,
    ) -> tuple[int, int, int]:
        """
        Parse InterProScan TSV output and enrich protein features.

        InterProScan TSV columns (0-indexed):
          0: protein accession (locus_tag)
          1: sequence MD5 digest
          2: sequence length
          3: analysis (e.g. Pfam, TIGRFAM, Gene3D)
          4: signature accession
          5: signature description
          6: start location
          7: stop location
          8: score (e-value)
          9: status (T = match)
         10: date
         11: IPR accession (or -)
         12: IPR description (or -)
         13: GO annotations (pipe-separated, or -)
         14: Pathways (or -)

        Returns (go_count, ipr_count, annotated_count).
        """
        go_count = 0
        ipr_count = 0
        annotated_tags: set[str] = set()

        if not tsv_path.exists() or tsv_path.stat().st_size == 0:
            return 0, 0, 0

        try:
            with open(tsv_path, newline="") as fh:
                reader = csv.reader(fh, delimiter="\t")
                for row in reader:
                    if len(row) < 12:
                        continue

                    locus_tag = row[0].strip()
                    analysis = row[3].strip()
                    sig_acc = row[4].strip()
                    sig_desc = row[5].strip()
                    ipr_acc = row[11].strip() if len(row) > 11 else "-"
                    ipr_desc = row[12].strip() if len(row) > 12 else "-"
                    go_str = row[13].strip() if len(row) > 13 else "-"
                    pathways = row[14].strip() if len(row) > 14 else "-"

                    feat = proteins.get(locus_tag)
                    if feat is None:
                        continue

                    annotated_tags.add(locus_tag)

                    # Add db_xref for the signature
                    db_xref = f"{analysis}:{sig_acc}"
                    if db_xref not in feat.db_xref:
                        feat.db_xref.append(db_xref)

                    # Add IPR ID
                    if ipr_acc and ipr_acc != "-":
                        if ipr_acc not in feat.ipr_ids:
                            feat.ipr_ids.append(ipr_acc)
                            ipr_count += 1

                        # If protein is still hypothetical, upgrade with IPR description
                        if feat.product == "hypothetical protein" and ipr_desc != "-":
                            feat.product = ipr_desc
                            feat.inference = f"protein motif:InterPro:{ipr_acc}"

                    # Add GO terms
                    if go_str and go_str != "-":
                        for go_term in go_str.split("|"):
                            go_term = go_term.strip()
                            if go_term.startswith("GO:") and go_term not in feat.go_terms:
                                feat.go_terms.append(go_term)
                                go_count += 1

                    # Add pathway to note
                    if pathways and pathways != "-":
                        pathway_note = f"pathway: {pathways}"
                        if feat.note:
                            if pathway_note not in feat.note:
                                feat.note += f"; {pathway_note}"
                        else:
                            feat.note = pathway_note

        except (csv.Error, OSError, ValueError) as e:
            logger.debug(f"Error parsing InterProScan output: {e}")

        return go_count, ipr_count, len(annotated_tags)
