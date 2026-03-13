"""
Checkpoint — Fossilization and resurrection of genome state.

Like amber preserving ancient organisms — checkpoints capture
the genome state at key moments so annotation can resume
after interruption without repeating expensive work.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

from darwin.rocks.models import (
    Contig,
    Feature,
    FeatureType,
    Genome,
    Strand,
)

logger = logging.getLogger("darwin.jar.checkpoint")


def save_checkpoint(
    genome: Genome,
    stage: str,
    output_dir: Path,
    config: dict | None = None,
    metadata: dict | None = None,
) -> Path:
    """
    Save genome state to a checkpoint file.

    The checkpoint captures all contigs, features, and annotations
    so work can resume from this point.
    """
    checkpoint_dir = Path(output_dir) / ".darwin_checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    checkpoint_path = checkpoint_dir / f"{stage}.json"

    data = {
        "stage": stage,
        "saved_at": time.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "genome": _serialize_genome(genome),
        "config": config or {},
        "metadata": metadata or {},
    }

    with open(checkpoint_path, "w") as fh:
        json.dump(data, fh, indent=2, default=str)

    logger.info(f"💾 Checkpoint saved: {stage} → {checkpoint_path}")
    return checkpoint_path


def load_checkpoint(output_dir: Path, stage: str | None = None) -> dict | None:
    """
    Load the latest (or specified) checkpoint.

    If stage is None, loads the most advanced checkpoint available.
    Returns None if no checkpoint exists.
    """
    checkpoint_dir = Path(output_dir) / ".darwin_checkpoints"
    if not checkpoint_dir.exists():
        return None

    if stage:
        path = checkpoint_dir / f"{stage}.json"
        if path.exists():
            return _load_file(path)
        return None

    # Find the most advanced checkpoint
    stage_order = [
        "genes_called",
        "proteins_found",
        "rrna_detected",
        "trna_detected",
        "crispr_detected",
        "signal_peptides_found",
        "operons_grouped",
        "taxonomy_inferred",
        "qc_completed",
        "annotation_ready",
    ]

    for s in reversed(stage_order):
        path = checkpoint_dir / f"{s}.json"
        if path.exists():
            return _load_file(path)

    return None


def list_checkpoints(output_dir: Path) -> list[dict]:
    """List all available checkpoints with metadata."""
    checkpoint_dir = Path(output_dir) / ".darwin_checkpoints"
    if not checkpoint_dir.exists():
        return []

    checkpoints = []
    for path in sorted(checkpoint_dir.glob("*.json")):
        try:
            with open(path) as fh:
                data = json.load(fh)
            checkpoints.append({
                "stage": data.get("stage", path.stem),
                "saved_at": data.get("saved_at", ""),
                "path": str(path),
            })
        except (json.JSONDecodeError, OSError):
            continue

    return checkpoints


def clear_checkpoints(output_dir: Path) -> int:
    """Remove all checkpoints. Returns number removed."""
    checkpoint_dir = Path(output_dir) / ".darwin_checkpoints"
    if not checkpoint_dir.exists():
        return 0

    count = 0
    for path in checkpoint_dir.glob("*.json"):
        path.unlink()
        count += 1

    if count:
        logger.info(f"🧹 Cleared {count} checkpoint(s)")
    return count


def _load_file(path: Path) -> dict:
    """Load and deserialize a checkpoint file."""
    with open(path) as fh:
        data: dict = json.load(fh)

    # Reconstruct genome from serialized data
    if "genome" in data:
        data["genome"] = _deserialize_genome(data["genome"])

    logger.info(f"💾 Checkpoint loaded: {data.get('stage', '?')} from {path}")
    return data


def _serialize_genome(genome: Genome) -> dict:
    """Serialize a Genome to a JSON-safe dict."""
    return {
        "name": genome.name,
        "organism": genome.organism,
        "taxonomy": genome.taxonomy,
        "source_file": str(genome.source_file) if genome.source_file else None,
        "contigs": [_serialize_contig(c) for c in genome.contigs],
    }


def _serialize_contig(contig: Contig) -> dict:
    return {
        "id": contig.id,
        "sequence": contig.sequence,
        "description": contig.description,
        "features": [_serialize_feature(f) for f in contig.features],
    }


def _serialize_feature(feature: Feature) -> dict:
    return {
        "type": feature.type.value,
        "start": feature.start,
        "end": feature.end,
        "strand": feature.strand.value,
        "score": feature.score,
        "contig_id": feature.contig_id,
        "locus_tag": feature.locus_tag,
        "product": feature.product,
        "inference": feature.inference,
        "translation": feature.translation,
        "gene": feature.gene,
        "note": feature.note,
        "db_xref": feature.db_xref,
    }


def _deserialize_genome(data: dict) -> Genome:
    """Reconstruct a Genome from serialized data."""
    contigs = [_deserialize_contig(c) for c in data.get("contigs", [])]
    return Genome(
        name=data["name"],
        contigs=contigs,
        source_file=Path(data["source_file"]) if data.get("source_file") else None,
        organism=data.get("organism", ""),
        taxonomy=data.get("taxonomy", ""),
    )


def _deserialize_contig(data: dict) -> Contig:
    features = [_deserialize_feature(f) for f in data.get("features", [])]
    return Contig(
        id=data["id"],
        sequence=data["sequence"],
        description=data.get("description", ""),
        features=features,
    )


def _deserialize_feature(data: dict) -> Feature:
    return Feature(
        type=FeatureType(data["type"]),
        start=data["start"],
        end=data["end"],
        strand=Strand(data.get("strand", ".")),
        score=data.get("score", 0.0),
        contig_id=data.get("contig_id", ""),
        locus_tag=data.get("locus_tag", ""),
        product=data.get("product", "hypothetical protein"),
        inference=data.get("inference", ""),
        translation=data.get("translation", ""),
        gene=data.get("gene", ""),
        note=data.get("note", ""),
        db_xref=data.get("db_xref", []),
    )
