"""
PhiSpyPlant — Prophage detection and annotation.

Feeds on: genome.loaded (raw sequence data)
Produces: prophages.detected (prophage region annotations)

Like bacteriophages lying dormant in a jar — integrated into the
host genome as prophages, waiting for the signal to wake up.
PhiSpy uses a combination of features (gene length, strand switching,
AT/GC skew, protein similarity) to detect these integrated viral
genomes lurking within the bacterial chromosome.

Uses PhiSpy for prophage region detection.
Install: conda install -c bioconda phispy
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

logger = logging.getLogger("darwin.flora.phispy")


class PhiSpyPlant(Organism):
    """PhiSpy prophage detection — finding dormant viruses."""

    name = "phispy"
    feeds_on_nutrients = [NutrientType.GENOME_LOADED]
    produces_nutrients = [NutrientType.PROPHAGES_DETECTED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        """Check if PhiSpy is available in soil."""
        return self.soil.has_phispy

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """
        Detect prophage regions in the genome.

        PhiSpy needs GenBank format input, so we build a minimal
        GenBank file from the genome FASTA, run PhiSpy, then parse
        the prophage coordinates from the output.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        cpus = config.get("cpus", 1)
        locus_prefix = config.get("locus_tag_prefix", "DARWIN")

        if not self.can_grow():
            self.logger.warning("🏜️ PhiSpy not in soil — cannot detect prophages")
            return Nutrient(
                type=NutrientType.PROPHAGES_DETECTED,
                data={
                    "genome": genome,
                    "prophage_count": 0,
                    "detection_skipped": True,
                    "config": config,
                },
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        self.logger.info(f"🦠 Scanning {genome.name} for prophage regions...")

        with tempfile.TemporaryDirectory(prefix="darwin_phispy_") as tmpdir:
            tmp = Path(tmpdir)
            input_gbk = tmp / "input.gbk"
            output_dir = tmp / "phispy_out"

            # Build minimal GenBank for PhiSpy
            self._write_minimal_genbank(genome, input_gbk)

            # Run PhiSpy
            cmd = [
                "PhiSpy.py",
                str(input_gbk),
                "-o", str(output_dir),
                "--threads", str(cpus),
                "--output_choice", "1",  # just coordinates
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
                    f"⚠️ PhiSpy failed (exit {proc.returncode}): {err_msg}"
                )
                return Nutrient(
                    type=NutrientType.PROPHAGES_DETECTED,
                    data={
                        "genome": genome,
                        "prophage_count": 0,
                        "detection_failed": True,
                        "config": config,
                    },
                    source=self.name,
                    correlation_id=nutrient.correlation_id,
                )

            prophage_count, status_counts = self._parse_phispy_output(
                genome, output_dir, locus_prefix
            )

        self.logger.info(
            f"🦠 Found {prophage_count} prophage region(s)"
            + (f" ({status_counts})" if status_counts else "")
        )

        return Nutrient(
            type=NutrientType.PROPHAGES_DETECTED,
            data={
                "genome": genome,
                "prophage_count": prophage_count,
                "status_counts": status_counts,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    @staticmethod
    def _write_minimal_genbank(genome: Genome, path: Path) -> None:
        """
        Write a minimal GenBank file from genome data.

        PhiSpy needs GenBank format. We create the simplest valid
        GenBank with just LOCUS, source feature, and sequence.
        """
        with open(path, "w") as fh:
            for contig in genome.contigs:
                fh.write(
                    f"LOCUS       {contig.id:<16} {contig.length} bp    "
                    f"DNA     linear   BCT\n"
                )
                fh.write(f"DEFINITION  {genome.name} {contig.description}\n")
                fh.write(f"ACCESSION   {contig.id}\n")
                fh.write(f"VERSION     {contig.id}\n")
                fh.write("FEATURES             Location/Qualifiers\n")
                fh.write(f"     source          1..{contig.length}\n")
                fh.write(f'                     /organism="{genome.organism or "Unknown"}"\n')
                fh.write('                     /mol_type="genomic DNA"\n')

                # Include CDS features if already called (helps PhiSpy accuracy)
                cds_features = contig.features_of_type(FeatureType.CDS)
                for f in sorted(cds_features, key=lambda x: x.start):
                    loc = f.location_str
                    fh.write(f"     CDS             {loc}\n")
                    if f.locus_tag:
                        fh.write(f'                     /locus_tag="{f.locus_tag}"\n')
                    fh.write(f'                     /product="{f.product}"\n')
                    if f.translation:
                        fh.write('                     /translation="')
                        seq = f.translation
                        fh.write(seq[:58])
                        for i in range(58, len(seq), 58):
                            fh.write(f"\n                     {seq[i:i+58]}")
                        fh.write('"\n')

                # Sequence
                fh.write("ORIGIN\n")
                seq = contig.sequence.lower()
                for i in range(0, len(seq), 60):
                    chunk = seq[i:i + 60]
                    parts = [chunk[j:j + 10] for j in range(0, len(chunk), 10)]
                    fh.write(f"{i + 1:>9} {' '.join(parts)}\n")
                fh.write("//\n")

    @staticmethod
    def _parse_phispy_output(
        genome: Genome,
        output_dir: Path,
        locus_prefix: str,
    ) -> tuple[int, dict[str, int]]:
        """
        Parse PhiSpy output and create PROPHAGE features.

        PhiSpy outputs prophage_coordinates.tsv with columns:
        Prophage_number, Contig, Start, Stop, ...

        Returns (prophage_count, status_counts_dict).
        """
        contig_map = {c.id: c for c in genome.contigs}
        prophage_count = 0
        status_counts: dict[str, int] = {}

        # PhiSpy output file — try multiple naming patterns
        candidates = [
            output_dir / "prophage_coordinates.tsv",
            output_dir / "prophage.tsv",
        ]
        # Also search for any TSV in output
        if output_dir.exists():
            candidates.extend(output_dir.glob("*prophage*"))
            candidates.extend(output_dir.glob("*.tsv"))

        coords_file = None
        for candidate in candidates:
            if isinstance(candidate, Path) and candidate.exists() and candidate.stat().st_size > 0:
                coords_file = candidate
                break

        if not coords_file:
            return 0, status_counts

        try:
            with open(coords_file, newline="") as fh:
                # PhiSpy TSV may or may not have headers
                first_line = fh.readline().strip()
                fh.seek(0)

                if first_line.startswith("Prophage") or first_line.startswith("#"):
                    reader = csv.DictReader(fh, delimiter="\t")
                    rows = list(reader)
                else:
                    # No header — use positional columns
                    # Format: pp_num, contig, start, stop, ...
                    rows = []
                    for line in fh:
                        parts = line.strip().split("\t")
                        if len(parts) >= 4:
                            rows.append({
                                "Prophage_number": parts[0],
                                "Contig": parts[1],
                                "Start": parts[2],
                                "Stop": parts[3],
                            })

                for row in rows:
                    # Try various column name patterns PhiSpy uses
                    contig_id = (
                        row.get("Contig", "")
                        or row.get("contig", "")
                        or row.get("Contig_id", "")
                    ).strip()

                    if not contig_id or contig_id not in contig_map:
                        # PhiSpy may use a different ID format; try matching
                        for cid in contig_map:
                            if cid in contig_id or contig_id in cid:
                                contig_id = cid
                                break
                        else:
                            continue

                    try:
                        start = int(row.get("Start", row.get("start", "0")))
                        end = int(row.get("Stop", row.get("stop", row.get("End", "0"))))
                    except ValueError:
                        continue

                    if start == 0 or end == 0:
                        continue

                    prophage_count += 1
                    locus_tag = f"{locus_prefix}_pp{prophage_count:04d}"

                    # Determine status/completeness if available
                    status = row.get("Status", row.get("status", "predicted")).strip()
                    if not status or status == "-":
                        status = "predicted"
                    status_counts[status] = status_counts.get(status, 0) + 1

                    length = abs(end - start)
                    note_parts = [
                        f"status: {status}",
                        f"length: {length}bp",
                    ]
                    pp_num = row.get("Prophage_number", row.get("pp_num", ""))
                    if pp_num:
                        note_parts.append(f"prophage ID: {pp_num}")

                    feature = Feature(
                        type=FeatureType.PROPHAGE,
                        start=min(start, end),
                        end=max(start, end),
                        strand=Strand.UNKNOWN,
                        contig_id=contig_id,
                        locus_tag=locus_tag,
                        product=f"{status} prophage region",
                        inference="ab initio prediction:PhiSpy",
                        note="; ".join(note_parts),
                    )
                    contig_map[contig_id].features.append(feature)

        except (csv.Error, KeyError, ValueError, OSError) as e:
            logger.debug(f"Error parsing PhiSpy output: {e}")

        return prophage_count, status_counts
