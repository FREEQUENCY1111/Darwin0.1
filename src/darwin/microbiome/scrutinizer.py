"""
Scrutinizer — Quality control decomposer.

Feeds on: ALL annotation signals (genes, tRNA, rRNA, proteins)
Produces: qc.completed

Like the cleanup crew in a jar ecosystem — it monitors the
health of the whole system. If something is off (too few genes,
missing tRNAs, weird density), it releases toxins (warnings)
into the water.

The Scrutinizer doesn't stop anything. It just reports.
Downstream organisms decide how to handle warnings.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

from darwin.flora.base import Organism
from darwin.rocks.models import FeatureType, Genome
from darwin.soil.nutrients import NutrientStore
from darwin.water.stream import Nutrient, NutrientType, Stream

logger = logging.getLogger("darwin.microbiome.scrutinizer")


@dataclass
class QCCheck:
    """Result of a single quality check."""

    name: str
    passed: bool
    value: float
    expected_range: str
    message: str
    severity: str = "info"  # info, warning, critical


class Scrutinizer(Organism):
    """
    Quality control — the ecosystem health monitor.

    Waits until it has enough signals to assess ecosystem health,
    then runs all QC checks and reports.
    """

    name = "scrutinizer"
    feeds_on_nutrients = [
        NutrientType.GENES_CALLED,
        NutrientType.PROTEINS_FOUND,
        NutrientType.TRNA_DETECTED,
        NutrientType.RRNA_DETECTED,
    ]
    produces_nutrients = [NutrientType.QC_COMPLETED]

    def __init__(self, stream: Stream, soil: NutrientStore) -> None:
        super().__init__(stream, soil)
        self._signals_received: dict[str, Nutrient] = {}
        self._expected_signals = {
            NutrientType.GENES_CALLED,
            NutrientType.TRNA_DETECTED,
            NutrientType.RRNA_DETECTED,
        }

    def plant(self) -> None:
        """Custom planting — we need to collect multiple signals."""
        for nt in self.feeds_on_nutrients:
            self.stream.subscribe(nt, self._collect_signal)
            self.logger.debug(f"🔬 Scrutinizer monitoring {nt.value}")

    async def _collect_signal(self, nutrient: Nutrient) -> None:
        """Collect signals and run QC when we have enough."""
        self._signals_received[nutrient.type.value] = nutrient
        self.logger.debug(
            f"🔬 Received {nutrient.type.value} "
            f"({len(self._signals_received)}/{len(self._expected_signals)} signals)"
        )

        # Check if we have enough signals to assess
        received_types = {NutrientType(k) for k in self._signals_received}
        if self._expected_signals.issubset(received_types):
            await self._run_qc(nutrient.correlation_id)

    async def _run_qc(self, correlation_id: str | None) -> None:
        """Run all quality checks on the assembled data."""
        # Get the genome from any signal (they all carry it)
        genome: Genome = list(self._signals_received.values())[0].data["genome"]

        self.logger.info(f"🔬 Running ecosystem health check on {genome.name}...")

        checks: list[QCCheck] = []

        # 1. Gene density check (~1 gene per kb)
        checks.append(self._check_gene_density(genome))

        # 2. tRNA complement check (expect ~20 amino acid types)
        checks.append(self._check_trna_complement(genome))

        # 3. rRNA operon check (5S + 16S + 23S)
        checks.append(self._check_rrna_operons(genome))

        # 4. Coding density check (~85-90% for prokaryotes)
        checks.append(self._check_coding_density(genome))

        # 5. Average gene length check (~950 bp)
        checks.append(self._check_avg_gene_length(genome))

        # 6. Overlapping features check
        checks.append(self._check_overlaps(genome))

        passed = sum(1 for c in checks if c.passed)
        total = len(checks)
        critical = [c for c in checks if c.severity == "critical" and not c.passed]

        self.logger.info(f"🔬 Health check: {passed}/{total} passed")
        for check in checks:
            icon = "✅" if check.passed else ("🔴" if check.severity == "critical" else "⚠️")
            self.logger.info(f"  {icon} {check.name}: {check.message}")

        # Release QC results into the water
        await self.stream.release(
            Nutrient(
                type=NutrientType.QC_COMPLETED,
                data={
                    "genome": genome,
                    "checks": [
                        {
                            "name": c.name,
                            "passed": c.passed,
                            "value": c.value,
                            "expected": c.expected_range,
                            "message": c.message,
                            "severity": c.severity,
                        }
                        for c in checks
                    ],
                    "passed": passed,
                    "total": total,
                    "critical_failures": len(critical),
                    "healthy": len(critical) == 0,
                    "config": list(self._signals_received.values())[0].data.get("config", {}),
                },
                source=self.name,
                correlation_id=correlation_id,
            )
        )

        # Release warnings for failed checks
        for check in checks:
            if not check.passed:
                await self.stream.release(
                    Nutrient(
                        type=NutrientType.WARNING,
                        data={"check": check.name, "message": check.message},
                        source=self.name,
                        correlation_id=correlation_id,
                    )
                )

    async def grow(self, nutrient: Nutrient) -> Nutrient | None:
        """Not used directly — we use _collect_signal instead."""
        return None

    def _check_gene_density(self, genome: Genome) -> QCCheck:
        cds = genome.features_of_type(FeatureType.CDS)
        density = len(cds) / (genome.total_length / 1000) if genome.total_length else 0
        passed = 0.5 <= density <= 1.5
        return QCCheck(
            name="gene_density",
            passed=passed,
            value=round(density, 3),
            expected_range="0.5-1.5 genes/kb",
            message=f"{density:.3f} genes/kb ({len(cds)} genes in {genome.total_length / 1000:.0f} kb)",
            severity="warning" if not passed else "info",
        )

    def _check_trna_complement(self, genome: Genome) -> QCCheck:
        trnas = genome.features_of_type(FeatureType.TRNA)
        amino_acids = set()
        for t in trnas:
            if "tRNA-" in t.product:
                aa = t.product.split("tRNA-")[1].split("(")[0]
                amino_acids.add(aa)
        count = len(amino_acids)
        passed = count >= 18  # expect ~20, but 18 is acceptable
        return QCCheck(
            name="trna_complement",
            passed=passed,
            value=count,
            expected_range=">=18 amino acid types",
            message=f"{count}/20 amino acid tRNAs found ({len(trnas)} total tRNAs)",
            severity="warning" if not passed else "info",
        )

    def _check_rrna_operons(self, genome: Genome) -> QCCheck:
        rrnas = genome.features_of_type(FeatureType.RRNA)
        types = {f.product for f in rrnas}
        has_5s = any("5S" in t for t in types)
        has_16s = any("16S" in t for t in types)
        has_23s = any("23S" in t for t in types)
        complete = has_5s and has_16s and has_23s
        found = [x for x in ["5S", "16S", "23S"] if any(x in t for t in types)]
        return QCCheck(
            name="rrna_operons",
            passed=complete,
            value=len(found),
            expected_range="5S + 16S + 23S",
            message=f"Found: {', '.join(found) if found else 'none'} ({len(rrnas)} total rRNAs)",
            severity="warning" if not complete else "info",
        )

    def _check_coding_density(self, genome: Genome) -> QCCheck:
        cds = genome.features_of_type(FeatureType.CDS)
        coding_bp = sum(f.length for f in cds)
        density = (coding_bp / genome.total_length * 100) if genome.total_length else 0
        passed = 70 <= density <= 95
        return QCCheck(
            name="coding_density",
            passed=passed,
            value=round(density, 1),
            expected_range="70-95%",
            message=f"{density:.1f}% coding ({coding_bp:,} bp of {genome.total_length:,} bp)",
            severity="warning" if not passed else "info",
        )

    def _check_avg_gene_length(self, genome: Genome) -> QCCheck:
        cds = genome.features_of_type(FeatureType.CDS)
        if not cds:
            return QCCheck(
                name="avg_gene_length",
                passed=False,
                value=0,
                expected_range="800-1200 bp",
                message="No CDS found",
                severity="critical",
            )
        avg = sum(f.length for f in cds) / len(cds)
        passed = 600 <= avg <= 1500
        return QCCheck(
            name="avg_gene_length",
            passed=passed,
            value=round(avg, 0),
            expected_range="600-1500 bp",
            message=f"{avg:.0f} bp average across {len(cds)} genes",
            severity="warning" if not passed else "info",
        )

    def _check_overlaps(self, genome: Genome) -> QCCheck:
        total_overlaps = 0
        for contig in genome.contigs:
            sorted_features = sorted(contig.features, key=lambda f: f.start)
            for i in range(1, len(sorted_features)):
                if sorted_features[i].start < sorted_features[i - 1].end:
                    total_overlaps += 1

        total_features = len(genome.all_features)
        overlap_pct = (total_overlaps / total_features * 100) if total_features else 0
        passed = overlap_pct < 15  # some overlap is normal
        return QCCheck(
            name="feature_overlaps",
            passed=passed,
            value=total_overlaps,
            expected_range="<15% of features",
            message=f"{total_overlaps} overlapping features ({overlap_pct:.1f}%)",
            severity="warning" if not passed else "info",
        )
