"""
PyhmmerPlant — Protein annotation organism.

Feeds on: genes.called (predicted CDS with translations)
Produces: proteins.found (annotated proteins with functions)

Like an epiphyte that grows on other plants — it takes the
raw gene predictions and enriches them with functional
annotation from HMM databases.

v0.2: Multi-database search with consensus — loops over ALL
HMM databases, tracks best hit per protein (lowest E-value wins),
stores db_xref cross-references.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from darwin.flora.base import Organism
from darwin.rocks.models import FeatureType, Genome
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.flora.pyhmmer")


@dataclass
class HMMHit:
    """Best HMM hit for a protein across all databases."""

    product: str
    db_name: str
    accession: str
    evalue: float
    score: float
    db_xref: str  # e.g. "Pfam:PF00001" or "TIGRFAMs:TIGR00123"


class PyhmmerPlant(Organism):
    """pyhmmer protein annotator — the epiphyte."""

    name = "pyhmmer"
    feeds_on_nutrients = [NutrientType.GENES_CALLED]
    produces_nutrients = [NutrientType.PROTEINS_FOUND]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)

    def can_grow(self) -> bool:
        return self.soil.has_hmm

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """
        Search protein translations against ALL HMM databases.

        Multi-database consensus: tracks best hit per protein
        across all databases (lowest E-value wins). Stores
        db_xref cross-references for each hit.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        evalue = config.get("evalue_threshold", 1e-10)
        cpus = config.get("cpus", 1)

        if not self.can_grow():
            self.logger.info("🏜️ No HMM databases in soil — passing through")
            return Nutrient(
                type=NutrientType.PROTEINS_FOUND,
                data={
                    "genome": genome,
                    "annotated_count": 0,
                    "total_cds": len(genome.features_of_type(FeatureType.CDS)),
                    "config": config,
                },
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        self.logger.info("🌿 Searching proteins against HMM databases...")

        cds_features = genome.features_of_type(FeatureType.CDS)
        proteins = {f.locus_tag: f.translation for f in cds_features if f.translation}

        if not proteins:
            self.logger.warning("No proteins to search — no translations available")
            return Nutrient(
                type=NutrientType.PROTEINS_FOUND,
                data={"genome": genome, "annotated_count": 0, "total_cds": 0, "config": config},
                source=self.name,
                correlation_id=nutrient.correlation_id,
            )

        # Track best hit per protein across ALL databases
        best_hits: dict[str, HMMHit] = {}

        try:
            import pyhmmer

            # Build digital sequences from translations
            alphabet = pyhmmer.easel.Alphabet.amino()
            sequences = []
            tag_map: dict[str | bytes, str] = {}

            for tag, seq in proteins.items():
                try:
                    ds = pyhmmer.easel.TextSequence(
                        name=tag.encode(),
                        sequence=seq,
                    ).digitize(alphabet)
                    sequences.append(ds)
                    # pyhmmer ≥0.12 returns str from hit.name;
                    # older versions return bytes — store both for safety
                    tag_map[tag] = tag
                    tag_map[tag.encode()] = tag
                except Exception:
                    continue

            if not sequences:
                self.logger.warning("No valid protein sequences to search")
                return Nutrient(
                    type=NutrientType.PROTEINS_FOUND,
                    data={
                        "genome": genome,
                        "annotated_count": 0,
                        "total_cds": len(cds_features),
                        "config": config,
                    },
                    source=self.name,
                    correlation_id=nutrient.correlation_id,
                )

            # Search ALL HMM databases and collect best hits
            hmm_databases = self.soil.get_hmm_databases()
            for hmm_db in hmm_databases:
                self.logger.info(f"  Searching {hmm_db.name}...")
                db_hits = 0

                with pyhmmer.plan7.HMMFile(str(hmm_db.path)) as hmm_file:
                    for hits in pyhmmer.hmmsearch(
                        hmm_file, sequences, E=evalue, cpus=cpus
                    ):  # type: ignore[arg-type]
                        for hit in hits:
                            if hit.included:
                                tag = tag_map.get(hit.name, "")  # type: ignore[call-overload]
                                if not tag:
                                    continue

                                # pyhmmer ≥0.12 returns str; older returns bytes
                                qname = hits.query.name  # type: ignore[union-attr]
                                product = qname.decode() if isinstance(qname, bytes) else qname

                                # Build accession and db_xref
                                qacc = getattr(hits.query, "accession", None)
                                if qacc:
                                    acc = qacc.decode() if isinstance(qacc, bytes) else qacc
                                else:
                                    acc = product

                                # Determine db_xref prefix from db name
                                db_prefix = hmm_db.name
                                if "pfam" in db_prefix.lower() or "pf" in db_prefix.lower():
                                    db_prefix = "Pfam"
                                elif "tigr" in db_prefix.lower():
                                    db_prefix = "TIGRFAMs"

                                db_xref = f"{db_prefix}:{acc}"

                                # Keep best hit per protein (lowest E-value)
                                existing = best_hits.get(tag)
                                if existing is None or hit.evalue < existing.evalue:
                                    best_hits[tag] = HMMHit(
                                        product=product,
                                        db_name=hmm_db.name,
                                        accession=acc,
                                        evalue=hit.evalue,
                                        score=hit.score,
                                        db_xref=db_xref,
                                    )
                                elif (
                                    existing is not None
                                    and hit.evalue <= existing.evalue * 10
                                    and db_xref != existing.db_xref
                                ):
                                    # Similar E-value from different DB — note both
                                    for f in cds_features:
                                        if f.locus_tag == tag:
                                            if db_xref not in f.db_xref:
                                                f.db_xref.append(db_xref)
                                            break

                                db_hits += 1

                self.logger.info(f"    {db_hits} included hits from {hmm_db.name}")

            # Apply best hits to features
            annotated_count = 0
            for tag, hmm_hit in best_hits.items():
                for f in cds_features:
                    if f.locus_tag == tag and f.product == "hypothetical protein":
                        f.product = hmm_hit.product
                        f.inference = f"protein motif:{hmm_hit.db_name}"
                        f.score = hmm_hit.score
                        if hmm_hit.db_xref not in f.db_xref:
                            f.db_xref.insert(0, hmm_hit.db_xref)
                        annotated_count += 1
                        break

        except ImportError:
            self.logger.warning("pyhmmer not installed — skipping HMM search")
            annotated_count = 0
        except Exception as e:
            self.logger.error(f"HMM search failed: {e}")
            annotated_count = 0

        hypothetical = sum(1 for f in cds_features if f.product == "hypothetical protein")
        self.logger.info(
            f"🧬 Annotated {annotated_count}/{len(cds_features)} proteins "
            f"({hypothetical} remain hypothetical)"
        )

        return Nutrient(
            type=NutrientType.PROTEINS_FOUND,
            data={
                "genome": genome,
                "annotated_count": annotated_count,
                "total_cds": len(cds_features),
                "hypothetical": hypothetical,
                "config": config,
            },
            source=self.name,
            correlation_id=nutrient.correlation_id,
        )
