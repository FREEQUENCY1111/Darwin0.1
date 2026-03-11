"""JSON output format writer — ideal for web API responses."""

from __future__ import annotations

import json
from pathlib import Path

from darwin.models import Genome


def write_json(genome: Genome, output: Path) -> Path:
    """Write genome annotations in a structured JSON format."""
    data = {
        "genome": genome.summary(),
        "contigs": [],
    }

    for contig in genome.contigs:
        contig_data = {
            "id": contig.id,
            "length": contig.length,
            "gc_content": round(contig.gc_content, 4),
            "features": [],
        }
        for feature in contig.features:
            feat_data = {
                "locus_tag": feature.locus_tag,
                "type": feature.feature_type.value,
                "start": feature.start,
                "end": feature.end,
                "strand": feature.strand.value,
                "length_bp": feature.length,
                "product": feature.product,
            }
            if feature.score is not None:
                feat_data["score"] = feature.score
            # Include extra attributes
            for key in ("inference", "anticodon", "evalue", "db_source"):
                if key in feature.attributes:
                    feat_data[key] = feature.attributes[key]

            contig_data["features"].append(feat_data)

        data["contigs"].append(contig_data)

    with open(output, "w") as fh:
        json.dump(data, fh, indent=2)

    return output


def genome_to_dict(genome: Genome) -> dict:
    """Convert genome to a JSON-serializable dict (for API responses)."""
    # Same structure as write_json but returns the dict
    data = {
        "genome": genome.summary(),
        "contigs": [],
    }
    for contig in genome.contigs:
        contig_data = {
            "id": contig.id,
            "length": contig.length,
            "gc_content": round(contig.gc_content, 4),
            "features": [
                {
                    "locus_tag": f.locus_tag,
                    "type": f.feature_type.value,
                    "start": f.start,
                    "end": f.end,
                    "strand": f.strand.value,
                    "length_bp": f.length,
                    "product": f.product,
                }
                for f in contig.features
            ],
        }
        data["contigs"].append(contig_data)
    return data
