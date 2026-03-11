"""
PyhmmerPlant — Protein annotation organism.

Feeds on: genes.called (predicted CDS with translations)
Produces: proteins.found (annotated proteins with functions)

Like an epiphyte that grows on other plants — it takes the
raw gene predictions and enriches them with functional
annotation from HMM databases.
"""

from __future__ import annotations

import logging

from darwin.flora.base import Organism
from darwin.rocks.models import FeatureType, Genome
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.flora.pyhmmer")


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
        Search protein translations against HMM databases.

        Takes CDS features with translations, runs hmmsearch,
        and updates product names for hits.
        """
        genome: Genome = nutrient.data["genome"]
        config = nutrient.data.get("config", {})
        evalue = config.get("evalue_threshold", 1e-10)

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

        annotated_count = 0

        try:
            import pyhmmer

            # Build digital sequences from translations
            alphabet = pyhmmer.easel.Alphabet.amino()
            sequences = []
            tag_map = {}

            for tag, seq in proteins.items():
                try:
                    ds = pyhmmer.easel.TextSequence(
                        name=tag.encode(),
                        sequence=seq,
                    ).digitize(alphabet)
                    sequences.append(ds)
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

            # Search each HMM database
            for hmm_db in self.soil.get_hmm_databases():
                self.logger.info(f"  Searching {hmm_db.name}...")

                with pyhmmer.plan7.HMMFile(str(hmm_db.path)) as hmm_file:
                    for hits in pyhmmer.hmmsearch(hmm_file, sequences, E=evalue):  # type: ignore[arg-type]
                        for hit in hits:
                            if hit.included:
                                tag = tag_map.get(hit.name, "")  # type: ignore[call-overload]
                                if tag:
                                    # Find and update the feature
                                    for f in cds_features:
                                        if (
                                            f.locus_tag == tag
                                            and f.product == "hypothetical protein"
                                        ):
                                            f.product = hits.query_name.decode()  # type: ignore[attr-defined]
                                            f.inference = f"protein motif:{hmm_db.name}"
                                            f.score = hit.score
                                            annotated_count += 1
                                            break

        except ImportError:
            self.logger.warning("pyhmmer not installed — skipping HMM search")
        except Exception as e:
            self.logger.error(f"HMM search failed: {e}")

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
