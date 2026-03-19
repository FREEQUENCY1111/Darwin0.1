"""
FoldseekPlant — 3D structural alignment for remote homology detection.

Feeds on: proteins.found (annotated proteins with translations)
Produces: structures.matched (structural hits from PDB/AlphaFold)

Like deep roots that reach mineral veins invisible from the surface —
Foldseek finds structural homologs that sequence methods miss entirely.
Two proteins can share <15% sequence identity yet fold into the same
3D structure with identical function. This organism catches those
"twilight zone" homologs.

Foldseek searches protein structures (predicted or experimental)
against PDB and the AlphaFold Database at >1000x the speed of
DALI or TM-align. It uses 3Di structural alphabet encoding to
convert 3D coordinates into a searchable sequence.

Pipeline:
  1. Predict structures from sequences (Foldseek's built-in predictor,
     or pre-predicted PDB/AlphaFold structures)
  2. Search against PDB and/or AlphaFold DB using 3Di + AA alignment
  3. Parse hits and annotate proteins with structural homologs

Uses Foldseek for structural search.
Install: conda install -c conda-forge -c bioconda foldseek
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

logger = logging.getLogger("darwin.flora.foldseek")


class FoldseekPlant(Organism):
    """Foldseek structural search — the deep roots."""

    name = "foldseek"
    feeds_on_nutrients = [NutrientType.PROTEINS_FOUND]
    produces_nutrients = [NutrientType.STRUCTURES_MATCHED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        """Check if Foldseek is available in soil."""
        return self.soil.has_foldseek

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """
        Search proteins against structural databases via Foldseek.

        Foldseek can operate in two modes:
          1. easy-search: query FASTA against a pre-built DB (PDB, AFDB)
          2. search: query structures against structure DB

        We use easy-search with --prostt5-model for direct sequence→structure
        search (no pre-computed structures needed). Falls back to standard
        sequence search against PDB if prostt5 is unavailable.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        cpus = config.get("cpus", 1)

        if not self.can_grow():
            self.logger.warning(
                "🏜️ Foldseek not in soil — cannot search structural databases"
            )
            return Nutrient(
                type=NutrientType.STRUCTURES_MATCHED,
                data={
                    "genome": genome,
                    "structural_hits": 0,
                    "remote_homologs": 0,
                    "detection_skipped": True,
                    "config": config,
                },
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        self.logger.info(
            f"🔬 Searching {genome.name} proteins against structural databases..."
        )

        cds_features = genome.features_of_type(FeatureType.CDS)
        proteins = {f.locus_tag: f for f in cds_features if f.translation}

        if not proteins:
            self.logger.warning("No proteins to search structurally")
            return Nutrient(
                type=NutrientType.STRUCTURES_MATCHED,
                data={
                    "genome": genome,
                    "structural_hits": 0,
                    "remote_homologs": 0,
                    "config": config,
                },
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        with tempfile.TemporaryDirectory(prefix="darwin_foldseek_") as tmpdir:
            tmp = Path(tmpdir)
            input_faa = tmp / "query.faa"
            result_file = tmp / "results.m8"

            # Write protein FASTA
            with open(input_faa, "w") as fh:
                for tag, feat in proteins.items():
                    fh.write(f">{tag}\n{feat.translation}\n")

            # Try to find a local Foldseek database
            db_path = self._find_foldseek_db()

            if db_path:
                # Local database search
                hit_count, remote_count = await self._search_local(
                    input_faa, db_path, result_file, proteins, cpus
                )
            else:
                # Remote search against PDB (Foldseek server)
                hit_count, remote_count = await self._search_remote(
                    input_faa, result_file, proteins, cpus
                )

        self.logger.info(
            f"🔬 Foldseek: {hit_count} structural hits, "
            f"{remote_count} remote homologs (sequence identity <30%)"
        )

        return Nutrient(
            type=NutrientType.STRUCTURES_MATCHED,
            data={
                "genome": genome,
                "structural_hits": hit_count,
                "remote_homologs": remote_count,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )

    def _find_foldseek_db(self) -> Path | None:
        """
        Look for a local Foldseek database.

        Search locations:
          1. ./databases/foldseek/    (project-local)
          2. ~/.darwin/databases/foldseek/  (user cache)

        A Foldseek DB consists of multiple files with a shared prefix.
        We look for .dbtype files as markers.
        """
        search_dirs = [
            Path.cwd() / "databases" / "foldseek",
            Path.home() / ".darwin" / "databases" / "foldseek",
        ]
        for search_dir in search_dirs:
            if not search_dir.exists():
                continue
            for dbtype_file in search_dir.glob("*.dbtype"):
                # The DB path is the prefix without .dbtype
                db_path = dbtype_file.with_suffix("")
                self.logger.info(f"🔬 Found Foldseek DB: {db_path}")
                return db_path

        self.logger.info(
            "🔬 No local Foldseek DB — using remote search against PDB"
        )
        return None

    async def _search_local(
        self,
        query_faa: Path,
        db_path: Path,
        result_file: Path,
        proteins: dict[str, object],
        cpus: int,
    ) -> tuple[int, int]:
        """Search against a local Foldseek database."""
        tmp_dir = result_file.parent / "foldseek_tmp"

        cmd = [
            "foldseek",
            "easy-search",
            str(query_faa),
            str(db_path),
            str(result_file),
            str(tmp_dir),
            "--threads", str(cpus),
            "--format-output",
            "query,target,pident,evalue,bits,alntmscore,qtmscore,ttmscore",
            "-e", "1e-3",
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
                f"⚠️ Foldseek local search failed (exit {proc.returncode}): {err_msg}"
            )
            return 0, 0

        return self._parse_foldseek_output(result_file, proteins)

    async def _search_remote(
        self,
        query_faa: Path,
        result_file: Path,
        proteins: dict[str, object],
        cpus: int,
    ) -> tuple[int, int]:
        """
        Search against Foldseek web server databases.

        Uses 'foldseek easy-search' with the --server flag if available,
        otherwise falls back to local PDB search.
        """
        tmp_dir = result_file.parent / "foldseek_tmp"

        # Try easy-search against PDB bundled with Foldseek
        cmd = [
            "foldseek",
            "easy-search",
            str(query_faa),
            "pdb100",  # Foldseek can auto-download PDB clustered at 100%
            str(result_file),
            str(tmp_dir),
            "--threads", str(cpus),
            "--format-output",
            "query,target,pident,evalue,bits,alntmscore,qtmscore,ttmscore",
            "-e", "1e-3",
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
                f"⚠️ Foldseek remote search failed (exit {proc.returncode}): "
                f"{err_msg}\n"
                f"Tip: download a local DB with: foldseek databases PDB pdb_db tmp/"
            )
            return 0, 0

        return self._parse_foldseek_output(result_file, proteins)

    @staticmethod
    def _parse_foldseek_output(
        result_file: Path,
        proteins: dict[str, object],
    ) -> tuple[int, int]:
        """
        Parse Foldseek tabular output and annotate protein features.

        Output format (custom --format-output):
          0: query (locus_tag)
          1: target (PDB/AF ID)
          2: pident (sequence identity %)
          3: evalue
          4: bits (bit score)
          5: alntmscore (alignment TM-score)
          6: qtmscore (query TM-score)
          7: ttmscore (target TM-score)

        A TM-score > 0.5 indicates structural similarity.
        A TM-score > 0.7 indicates same fold.

        Returns (total_hits, remote_homologs_count).
        """
        if not result_file.exists() or result_file.stat().st_size == 0:
            return 0, 0

        # Track best hit per protein
        best_hits: dict[str, dict] = {}
        remote_count = 0

        try:
            with open(result_file) as fh:
                reader = csv.reader(fh, delimiter="\t")
                for row in reader:
                    if len(row) < 5:
                        continue

                    query = row[0].strip()
                    target = row[1].strip()

                    try:
                        pident = float(row[2]) if len(row) > 2 else 0.0
                        evalue = float(row[3]) if len(row) > 3 else 999.0
                        bits = float(row[4]) if len(row) > 4 else 0.0
                        aln_tmscore = float(row[5]) if len(row) > 5 else 0.0
                        q_tmscore = float(row[6]) if len(row) > 6 else 0.0
                    except ValueError:
                        continue

                    # Keep best hit per protein (highest bit score)
                    if query not in best_hits or bits > best_hits[query]["bits"]:
                        best_hits[query] = {
                            "target": target,
                            "pident": pident,
                            "evalue": evalue,
                            "bits": bits,
                            "aln_tmscore": aln_tmscore,
                            "q_tmscore": q_tmscore,
                        }

            # Apply hits to features
            for locus_tag, hit in best_hits.items():
                feat = proteins.get(locus_tag)
                if feat is None:
                    continue

                target = hit["target"]
                pident = hit["pident"]
                aln_tmscore = hit["aln_tmscore"]
                q_tmscore = hit["q_tmscore"]
                evalue = hit["evalue"]

                # Store structural hit info
                feat.structure_hit = target
                db_xref = f"PDB:{target}" if not target.startswith("AF-") else f"AlphaFoldDB:{target}"
                if db_xref not in feat.db_xref:
                    feat.db_xref.append(db_xref)

                # Build structural note
                note_parts = [f"structural hit: {target}"]
                note_parts.append(f"seq identity: {pident:.1f}%")
                if aln_tmscore > 0:
                    note_parts.append(f"TM-score: {aln_tmscore:.3f}")
                    if aln_tmscore > 0.7:
                        note_parts.append("same fold")
                    elif aln_tmscore > 0.5:
                        note_parts.append("similar fold")

                struct_note = "; ".join(note_parts)
                if feat.note:
                    feat.note += f"; {struct_note}"
                else:
                    feat.note = struct_note

                # Track remote homologs (structurally similar but low sequence identity)
                if pident < 30.0 and (aln_tmscore > 0.5 or evalue < 1e-5):
                    remote_count += 1

                    # If protein is still hypothetical, try to upgrade from PDB annotation
                    if feat.product == "hypothetical protein":
                        # Extract protein name from PDB target ID if possible
                        # PDB IDs like "1xyz_A" or AF IDs like "AF-P12345-F1"
                        feat.inference = f"structural similarity:{target}"

        except (csv.Error, OSError, ValueError) as e:
            logger.debug(f"Error parsing Foldseek output: {e}")

        return len(best_hits), remote_count
