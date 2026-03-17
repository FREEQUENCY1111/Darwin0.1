"""
GeccoPlant — Biosynthetic gene cluster (BGC) detection.

Feeds on: genome.loaded (raw sequence data)
Produces: bgc.detected (BGC region annotations)

Like specialized microbes that produce secondary metabolites —
antibiotics, siderophores, toxins — the Gecco plant identifies
biosynthetic gene clusters that encode the molecular factories
for these natural products. Finding these clusters reveals the
organism's chemical arsenal and biotechnological potential.

Uses GECCO (GEne Cluster prediction with Conditional random fields;
Optimized) for ML-based BGC detection.
Install: conda install -c bioconda gecco
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

logger = logging.getLogger("darwin.flora.gecco")


class GeccoPlant(Organism):
    """GECCO BGC detection — finding natural product factories."""

    name = "gecco"
    feeds_on_nutrients = [NutrientType.GENOME_LOADED]
    produces_nutrients = [NutrientType.BGC_DETECTED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        """Check if GECCO is available in soil."""
        return self.soil.has_gecco

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """
        Detect biosynthetic gene clusters in the genome.

        Runs GECCO on the genome FASTA, parses the cluster TSV output,
        and creates BGC features for each predicted cluster with type
        classification and confidence scores.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        cpus = config.get("cpus", 1)
        locus_prefix = config.get("locus_tag_prefix", "DARWIN")

        if not self.can_grow():
            self.logger.warning("🏜️ GECCO not in soil — cannot detect BGCs")
            return Nutrient(
                type=NutrientType.BGC_DETECTED,
                data={
                    "genome": genome,
                    "bgc_count": 0,
                    "bgc_types": {},
                    "detection_skipped": True,
                    "config": config,
                },
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        self.logger.info(
            f"🧪 Scanning {genome.name} for biosynthetic gene clusters..."
        )

        with tempfile.TemporaryDirectory(prefix="darwin_gecco_") as tmpdir:
            tmp = Path(tmpdir)
            input_fasta = tmp / "input.fasta"
            output_dir = tmp / "gecco_out"

            # Write genome to temp FASTA
            with open(input_fasta, "w") as fh:
                for contig in genome.contigs:
                    fh.write(f">{contig.id}\n{contig.sequence}\n")

            # Run GECCO
            cmd = [
                "gecco",
                "run",
                "--genome", str(input_fasta),
                "--output-dir", str(output_dir),
                "--jobs", str(cpus),
            ]

            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, stderr = await proc.communicate()

            if proc.returncode != 0:
                err_msg = stderr.decode()[:500] or stdout.decode()[:500]
                # Negative exit codes mean killed by signal (e.g. -4 = SIGILL)
                if proc.returncode < 0:
                    import signal as _sig
                    try:
                        sig_name = _sig.Signals(-proc.returncode).name
                    except (ValueError, AttributeError):
                        sig_name = f"signal {-proc.returncode}"
                    self.logger.warning(
                        f"⚠️ GECCO crashed ({sig_name}) — this may be a "
                        f"platform compatibility issue (e.g. ARM Mac). "
                        f"Try: conda install -c conda-forge -c bioconda gecco"
                    )
                else:
                    self.logger.warning(
                        f"⚠️ GECCO failed (exit {proc.returncode}): {err_msg}"
                    )
                return Nutrient(
                    type=NutrientType.BGC_DETECTED,
                    data={
                        "genome": genome,
                        "bgc_count": 0,
                        "bgc_types": {},
                        "detection_failed": True,
                        "config": config,
                    },
                    source=self.name,
                    correlation_id=nutrient.correlation_id,
                )

            bgc_count, bgc_types = self._parse_gecco_output(
                genome, output_dir, locus_prefix
            )

        type_summary = ", ".join(
            f"{k}: {v}" for k, v in sorted(bgc_types.items())
        )
        self.logger.info(
            f"🧪 Found {bgc_count} biosynthetic gene cluster(s)"
            + (f" ({type_summary})" if type_summary else "")
        )

        return Nutrient(
            type=NutrientType.BGC_DETECTED,
            data={
                "genome": genome,
                "bgc_count": bgc_count,
                "bgc_types": bgc_types,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    @staticmethod
    def _parse_gecco_output(
        genome: Genome,
        output_dir: Path,
        locus_prefix: str,
    ) -> tuple[int, dict[str, int]]:
        """
        Parse GECCO output and create BGC features.

        GECCO outputs *.clusters.tsv with columns including:
        sequence_id, bgc_id, start, end, average_p, type, ...

        Also tries *.features.tsv as fallback.

        Returns (bgc_count, type_distribution_dict).
        """
        contig_map = {c.id: c for c in genome.contigs}
        bgc_count = 0
        bgc_types: dict[str, int] = {}

        if not output_dir.exists():
            return 0, bgc_types

        # Find cluster output file — GECCO uses various naming patterns
        cluster_files = list(output_dir.glob("*.clusters.tsv"))
        if not cluster_files:
            cluster_files = list(output_dir.glob("*cluster*"))
        if not cluster_files:
            cluster_files = list(output_dir.glob("*.tsv"))

        for cluster_file in cluster_files:
            if not cluster_file.exists() or cluster_file.stat().st_size == 0:
                continue

            try:
                with open(cluster_file, newline="") as fh:
                    reader = csv.DictReader(fh, delimiter="\t")
                    for row in reader:
                        # Try various column name patterns
                        seq_id = (
                            row.get("sequence_id", "")
                            or row.get("sequence", "")
                            or row.get("seqid", "")
                        ).strip()

                        if not seq_id:
                            continue

                        # GECCO may prefix contig IDs; try to match
                        matched_id = None
                        if seq_id in contig_map:
                            matched_id = seq_id
                        else:
                            for cid in contig_map:
                                if cid in seq_id or seq_id in cid:
                                    matched_id = cid
                                    break

                        if not matched_id:
                            continue

                        try:
                            start = int(row.get("start", row.get("Start", "0")))
                            end = int(row.get("end", row.get("End", "0")))
                        except ValueError:
                            continue

                        if start == 0 or end == 0:
                            continue

                        bgc_count += 1
                        locus_tag = f"{locus_prefix}_bgc{bgc_count:04d}"

                        # Cluster type
                        cluster_type = (
                            row.get("type", "")
                            or row.get("Type", "")
                            or row.get("product_type", "")
                            or "Unknown"
                        ).strip()

                        bgc_types[cluster_type] = bgc_types.get(cluster_type, 0) + 1

                        # Confidence score
                        avg_p = row.get("average_p", row.get("score", "")).strip()
                        confidence = ""
                        if avg_p:
                            try:
                                confidence = f"{float(avg_p):.3f}"
                            except ValueError:
                                confidence = avg_p

                        # Number of genes in cluster
                        n_genes = row.get("n_genes", row.get("genes", "")).strip()

                        # Build note
                        note_parts = [f"type: {cluster_type}"]
                        if confidence:
                            note_parts.append(f"confidence: {confidence}")
                        if n_genes:
                            note_parts.append(f"genes in cluster: {n_genes}")

                        bgc_id = row.get("bgc_id", row.get("cluster_id", "")).strip()
                        if bgc_id:
                            note_parts.append(f"cluster ID: {bgc_id}")

                        length = abs(end - start)
                        note_parts.append(f"length: {length}bp")

                        feature = Feature(
                            type=FeatureType.BGC,
                            start=min(start, end),
                            end=max(start, end),
                            strand=Strand.UNKNOWN,
                            contig_id=matched_id,
                            locus_tag=locus_tag,
                            product=f"{cluster_type} biosynthetic gene cluster",
                            inference="ab initio prediction:GECCO",
                            note="; ".join(note_parts),
                        )
                        contig_map[matched_id].features.append(feature)

            except (csv.Error, KeyError, ValueError, OSError) as e:
                logger.debug(f"Error parsing GECCO output {cluster_file.name}: {e}")
                continue

        return bgc_count, bgc_types
