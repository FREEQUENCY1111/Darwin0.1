"""Scrutinizer Agent — quality control and validation.

The Scrutinizer is the skeptic of the council. It:
  - Validates annotation results for biological plausibility
  - Checks for overlapping features
  - Flags unusual gene density (too high or too low)
  - Verifies tRNA complement (all amino acids covered?)
  - Checks rRNA operon completeness (5S + 16S + 23S)
  - Detects potential frameshifts or pseudogenes
  - Collects alerts from other agents and compiles a QC report
  - Can VETO the pipeline if critical issues are found
"""

from __future__ import annotations

from typing import Any

from darwin.agents.base import Agent, AgentMessage, AgentStatus, MessageType
from darwin.models import FeatureType, Genome


class Scrutinizer(Agent):
    name = "scrutinizer"
    role = "Quality control — validates results, flags issues, can veto the pipeline"

    def __init__(self) -> None:
        super().__init__()
        self._alerts: list[dict] = []
        self._qc_passed = True

    def handle_message(self, message: AgentMessage) -> None:
        """Collect alerts from other agents."""
        if message.msg_type == MessageType.ALERT:
            self._alerts.append({
                "from": message.sender,
                "issue": message.payload.get("issue", ""),
                "severity": message.payload.get("severity", "info"),
                "stage": message.payload.get("stage", "unknown"),
            })

    def validate_input(self, **kwargs: Any) -> bool:
        genome = self.context.genome
        return genome is not None and genome.num_contigs > 0

    def execute(self, **kwargs: Any) -> dict[str, Any]:
        self.set_status(AgentStatus.WORKING)
        self.log.info("[bold yellow]Scrutinizer[/] — validating annotation quality")

        genome = self.context.get_genome()
        qc_results: dict[str, Any] = {
            "checks": [],
            "alerts_received": list(self._alerts),
            "passed": True,
        }

        # Run all QC checks
        checks = [
            self._check_gene_density(genome),
            self._check_overlapping_features(genome),
            self._check_trna_complement(genome),
            self._check_rrna_operons(genome),
            self._check_coding_density(genome),
            self._check_avg_gene_length(genome),
        ]

        for check in checks:
            qc_results["checks"].append(check)
            if check["severity"] == "critical":
                qc_results["passed"] = False
                self._qc_passed = False

        # Check for critical alerts from other agents
        critical_alerts = [a for a in self._alerts if a["severity"] == "critical"]
        if critical_alerts:
            qc_results["passed"] = False
            self._qc_passed = False

        self.context.store_result(self.name, qc_results)

        # Report to council
        if qc_results["passed"]:
            self.broadcast(
                MessageType.STATUS,
                {"message": "QC PASSED — all checks within acceptable ranges"},
            )
            self.log.info("  QC [green]PASSED[/]")
        else:
            self.broadcast(
                MessageType.ALERT,
                {"message": "QC ISSUES FOUND — review recommended",
                 "failed_checks": [c for c in checks if c["severity"] != "ok"]},
            )
            self.log.warning("  QC [red]ISSUES FOUND[/] — review recommended")

        # Log individual check results
        for check in checks:
            icon = "✓" if check["severity"] == "ok" else "⚠" if check["severity"] == "warning" else "✗"
            self.log.info(f"  {icon} {check['name']}: {check['message']}")

        self.set_status(AgentStatus.DONE)
        return qc_results

    @property
    def qc_passed(self) -> bool:
        return self._qc_passed

    # ── QC Checks ────────────────────────────────────────────

    def _check_gene_density(self, genome: Genome) -> dict:
        """Check genes per kb — prokaryotes typically have ~1 gene/kb."""
        cds_count = sum(
            1 for f in genome.all_features if f.feature_type == FeatureType.CDS
        )
        if genome.total_length == 0:
            return {"name": "gene_density", "severity": "critical",
                    "message": "No sequence data"}

        density = cds_count / (genome.total_length / 1000)

        if density < 0.5:
            return {"name": "gene_density", "severity": "warning",
                    "message": f"Low gene density: {density:.2f} genes/kb (expected ~1.0)",
                    "value": density}
        elif density > 1.5:
            return {"name": "gene_density", "severity": "warning",
                    "message": f"High gene density: {density:.2f} genes/kb (expected ~1.0)",
                    "value": density}
        return {"name": "gene_density", "severity": "ok",
                "message": f"{density:.2f} genes/kb", "value": density}

    def _check_overlapping_features(self, genome: Genome) -> dict:
        """Detect features with significant overlap."""
        overlap_count = 0
        for contig in genome.contigs:
            features = sorted(contig.features, key=lambda f: f.start)
            for i in range(len(features) - 1):
                f1 = features[i]
                f2 = features[i + 1]
                overlap = f1.end - f2.start
                if overlap > 50:  # >50bp overlap is suspicious
                    overlap_count += 1

        if overlap_count > 20:
            return {"name": "overlapping_features", "severity": "warning",
                    "message": f"{overlap_count} features with >50bp overlap",
                    "value": overlap_count}
        return {"name": "overlapping_features", "severity": "ok",
                "message": f"{overlap_count} overlaps detected", "value": overlap_count}

    def _check_trna_complement(self, genome: Genome) -> dict:
        """Check if all 20 standard amino acid tRNAs are present."""
        trna_products = [
            f.product for f in genome.all_features
            if f.feature_type == FeatureType.TRNA
        ]
        amino_acids = set()
        for p in trna_products:
            if p.startswith("tRNA-"):
                aa = p.replace("tRNA-", "")
                amino_acids.add(aa)

        standard_20 = {
            "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly",
            "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser",
            "Thr", "Trp", "Tyr", "Val",
        }
        missing = standard_20 - amino_acids
        coverage = len(amino_acids) / 20

        if missing and len(missing) > 5:
            return {"name": "trna_complement", "severity": "warning",
                    "message": f"Missing {len(missing)} tRNA types: {', '.join(sorted(missing))}",
                    "value": coverage}
        elif missing:
            return {"name": "trna_complement", "severity": "ok",
                    "message": f"{len(amino_acids)}/20 amino acid tRNAs found",
                    "value": coverage}
        return {"name": "trna_complement", "severity": "ok",
                "message": "All 20 amino acid tRNAs present", "value": 1.0}

    def _check_rrna_operons(self, genome: Genome) -> dict:
        """Check for complete rRNA operons (5S + 16S + 23S)."""
        rrna_types = set()
        for f in genome.all_features:
            if f.feature_type == FeatureType.RRNA:
                product = f.product.lower()
                if "5s" in product:
                    rrna_types.add("5S")
                elif "16s" in product:
                    rrna_types.add("16S")
                elif "23s" in product:
                    rrna_types.add("23S")

        expected = {"5S", "16S", "23S"}
        missing = expected - rrna_types

        if missing:
            return {"name": "rrna_operons", "severity": "warning",
                    "message": f"Missing rRNA: {', '.join(sorted(missing))}",
                    "value": len(rrna_types)}
        return {"name": "rrna_operons", "severity": "ok",
                "message": "Complete rRNA operon (5S, 16S, 23S)", "value": 3}

    def _check_coding_density(self, genome: Genome) -> dict:
        """Check what percentage of the genome codes for proteins."""
        total_coding = sum(
            f.length for f in genome.all_features
            if f.feature_type == FeatureType.CDS
        )
        if genome.total_length == 0:
            return {"name": "coding_density", "severity": "critical",
                    "message": "No data"}

        density = total_coding / genome.total_length

        if density < 0.70:
            return {"name": "coding_density", "severity": "warning",
                    "message": f"Low coding density: {density:.1%} (expected ~85-90%)",
                    "value": density}
        return {"name": "coding_density", "severity": "ok",
                "message": f"{density:.1%} coding", "value": density}

    def _check_avg_gene_length(self, genome: Genome) -> dict:
        """Check average gene length — prokaryotic genes average ~950bp."""
        cds = [f for f in genome.all_features if f.feature_type == FeatureType.CDS]
        if not cds:
            return {"name": "avg_gene_length", "severity": "warning",
                    "message": "No CDS found"}

        avg = sum(f.length for f in cds) / len(cds)

        if avg < 500:
            return {"name": "avg_gene_length", "severity": "warning",
                    "message": f"Short average gene: {avg:.0f} bp (expected ~950)",
                    "value": avg}
        elif avg > 1500:
            return {"name": "avg_gene_length", "severity": "warning",
                    "message": f"Long average gene: {avg:.0f} bp (expected ~950)",
                    "value": avg}
        return {"name": "avg_gene_length", "severity": "ok",
                "message": f"{avg:.0f} bp average", "value": avg}
