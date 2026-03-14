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
from darwin.rocks.models import FeatureType, Genome, Strand
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
        NutrientType.PLASMIDS_CLASSIFIED,
        NutrientType.MOBILE_ELEMENTS_FOUND,
        NutrientType.RESISTANCE_GENES_FOUND,
        NutrientType.PROPHAGES_DETECTED,
        NutrientType.BGC_DETECTED,
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
        # Optionally wait for tool-dependent signals
        if self.soil.has_mob_suite:
            self._expected_signals.add(NutrientType.PLASMIDS_CLASSIFIED)
        if self.soil.has_isescan:
            self._expected_signals.add(NutrientType.MOBILE_ELEMENTS_FOUND)
        if self.soil.has_abricate:
            self._expected_signals.add(NutrientType.RESISTANCE_GENES_FOUND)
        if self.soil.has_phispy:
            self._expected_signals.add(NutrientType.PROPHAGES_DETECTED)
        if self.soil.has_gecco:
            self._expected_signals.add(NutrientType.BGC_DETECTED)

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

        # 7. rRNA length validation
        checks.append(self._check_rrna_lengths(genome))

        # 8. Start codon distribution
        checks.append(self._check_start_codon_distribution(genome))

        # 9. Strand bias
        checks.append(self._check_strand_bias(genome))

        # 10. rRNA copy number
        checks.append(self._check_rrna_copy_number(genome))

        # 11. IS element density (if ISEScan ran)
        if NutrientType.MOBILE_ELEMENTS_FOUND.value in self._signals_received:
            checks.append(self._check_is_element_density(genome))

        # 12. Plasmid classification report (if MOB-suite ran)
        if NutrientType.PLASMIDS_CLASSIFIED.value in self._signals_received:
            checks.append(self._check_plasmid_classification(genome))

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

    def _check_rrna_lengths(self, genome: Genome) -> QCCheck:
        """Validate rRNA feature lengths against expected ranges."""
        rrnas = genome.features_of_type(FeatureType.RRNA)
        if not rrnas:
            return QCCheck(
                name="rrna_lengths",
                passed=True,
                value=0,
                expected_range="16S:1400-1600, 23S:2800-3000, 5S:100-150",
                message="No rRNAs to validate",
                severity="info",
            )

        expected_ranges = {
            "16S": (1400, 1600),
            "23S": (2800, 3000),
            "5S": (100, 150),
        }
        issues = []
        checked = 0
        for f in rrnas:
            for rtype, (rmin, rmax) in expected_ranges.items():
                if rtype in f.product:
                    checked += 1
                    if not (rmin <= f.length <= rmax):
                        issues.append(f"{rtype}={f.length}bp (expect {rmin}-{rmax})")

        passed = len(issues) == 0
        return QCCheck(
            name="rrna_lengths",
            passed=passed,
            value=checked - len(issues),
            expected_range="16S:1400-1600, 23S:2800-3000, 5S:100-150",
            message=f"{checked} rRNAs checked, {len(issues)} out of range"
            + (f": {', '.join(issues)}" if issues else ""),
            severity="warning" if not passed else "info",
        )

    def _check_start_codon_distribution(self, genome: Genome) -> QCCheck:
        """Check ATG/GTG/TTG start codon ratios (expect ~80/10/10)."""
        cds = genome.features_of_type(FeatureType.CDS)
        if not cds:
            return QCCheck(
                name="start_codon_dist",
                passed=True,
                value=0,
                expected_range="ATG ~80%, GTG ~10%, TTG ~10%",
                message="No CDS to check",
                severity="info",
            )

        counts: dict[str, int] = {"ATG": 0, "GTG": 0, "TTG": 0, "other": 0}
        for f in cds:
            # Extract start codon from the contig sequence
            for contig in genome.contigs:
                if contig.id == f.contig_id:
                    if f.strand == Strand.FORWARD:
                        codon = contig.sequence[f.start - 1 : f.start + 2].upper()
                    else:
                        # Reverse complement of last 3 bases
                        nuc = contig.sequence[f.end - 3 : f.end]
                        comp = str.maketrans("ATCGatcg", "TAGCtagc")
                        codon = nuc.translate(comp)[::-1].upper()
                    if codon in counts:
                        counts[codon] += 1
                    else:
                        counts["other"] += 1
                    break

        total = sum(counts.values())
        if total == 0:
            return QCCheck(
                name="start_codon_dist",
                passed=True,
                value=0,
                expected_range="ATG ~80%, GTG ~10%, TTG ~10%",
                message="Could not extract start codons",
                severity="info",
            )

        atg_pct = counts["ATG"] / total * 100
        # ATG should dominate (~80%), but 50-95% is acceptable
        passed = 50 <= atg_pct <= 95
        dist = ", ".join(f"{k}:{v}" for k, v in counts.items() if v > 0)
        return QCCheck(
            name="start_codon_dist",
            passed=passed,
            value=round(atg_pct, 1),
            expected_range="ATG 50-95%",
            message=f"ATG={atg_pct:.1f}% ({dist})",
            severity="warning" if not passed else "info",
        )

    def _check_strand_bias(self, genome: Genome) -> QCCheck:
        """Check CDS should be roughly balanced across strands (40-60%)."""
        cds = genome.features_of_type(FeatureType.CDS)
        if not cds:
            return QCCheck(
                name="strand_bias",
                passed=True,
                value=0,
                expected_range="40-60% per strand",
                message="No CDS to check",
                severity="info",
            )

        fwd = sum(1 for f in cds if f.strand == Strand.FORWARD)
        fwd_pct = fwd / len(cds) * 100
        passed = 30 <= fwd_pct <= 70  # generous range
        return QCCheck(
            name="strand_bias",
            passed=passed,
            value=round(fwd_pct, 1),
            expected_range="30-70% forward strand",
            message=f"{fwd_pct:.1f}% forward ({fwd}/{len(cds)})",
            severity="warning" if not passed else "info",
        )

    def _check_rrna_copy_number(self, genome: Genome) -> QCCheck:
        """Check 16S rRNA copy number vs expected for genome size."""
        rrnas = genome.features_of_type(FeatureType.RRNA)
        copies_16s = sum(1 for f in rrnas if "16S" in f.product)

        # Expected copies based on genome size (rough guide)
        mb = genome.total_length / 1_000_000
        if mb < 1:
            expected = "1-2"
            passed = 1 <= copies_16s <= 3
        elif mb < 3:
            expected = "1-7"
            passed = 1 <= copies_16s <= 10
        elif mb < 6:
            expected = "3-10"
            passed = 1 <= copies_16s <= 15
        else:
            expected = "5-15"
            passed = 1 <= copies_16s <= 20

        return QCCheck(
            name="rrna_copy_number",
            passed=passed,
            value=copies_16s,
            expected_range=f"{expected} copies for {mb:.1f} Mb genome",
            message=f"{copies_16s} copies of 16S rRNA for {mb:.1f} Mb genome",
            severity="warning" if not passed else "info",
        )

    def _check_is_element_density(self, genome: Genome) -> QCCheck:
        """Check IS element density — warn if >50/Mb (possible misassembly)."""
        is_elements = genome.features_of_type(FeatureType.MOBILE_ELEMENT)
        mb = genome.total_length / 1_000_000
        density = len(is_elements) / mb if mb > 0 else 0
        passed = density <= 50
        return QCCheck(
            name="is_element_density",
            passed=passed,
            value=round(density, 1),
            expected_range="<=50 IS/Mb",
            message=f"{density:.1f} IS elements/Mb ({len(is_elements)} in {mb:.1f} Mb)",
            severity="warning" if not passed else "info",
        )

    def _check_plasmid_classification(self, genome: Genome) -> QCCheck:
        """Report plasmid classification results."""
        plasmids = [c for c in genome.contigs if c.replicon_type == "plasmid"]
        chromosomes = [c for c in genome.contigs if c.replicon_type == "chromosome"]
        unclassified = [c for c in genome.contigs if not c.replicon_type]
        total = genome.num_contigs
        # Always passes — informational check
        return QCCheck(
            name="plasmid_classification",
            passed=True,
            value=len(plasmids),
            expected_range="informational",
            message=(
                f"{len(plasmids)} plasmid(s), {len(chromosomes)} chromosomal, "
                f"{len(unclassified)} unclassified out of {total} contigs"
            ),
            severity="info",
        )
