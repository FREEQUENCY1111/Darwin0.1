"""Code-checker Chris — pipeline integrity and validation agent.

Chris ensures the pipeline itself is running correctly:
  - Validates that tool outputs are well-formed (valid GFF, FASTA, etc.)
  - Checks that intermediate files exist and are non-empty
  - Verifies locus tag uniqueness and consistency
  - Ensures output files are valid and parseable
  - Detects silent failures (tools that exit 0 but produce bad output)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from darwin.agents.base import Agent, AgentStatus, MessageType


class CodeCheckerChris(Agent):
    name = "code_checker_chris"
    role = "Validates pipeline integrity — checks outputs are well-formed and consistent"

    def validate_input(self, **kwargs: Any) -> bool:
        return self.context.genome is not None

    def execute(self, **kwargs: Any) -> dict[str, Any]:
        self.set_status(AgentStatus.WORKING)
        self.log.info("[bold red]Code-checker Chris[/] — validating pipeline integrity")

        genome = self.context.get_genome()
        checks: list[dict] = []

        # 1. Verify locus tag uniqueness
        checks.append(self._check_locus_tags(genome))

        # 2. Verify feature coordinates are valid
        checks.append(self._check_coordinates(genome))

        # 3. Verify output files exist and are non-empty
        synth_result = self.context.get_result("synthesizer")
        if synth_result:
            checks.append(self._check_output_files(synth_result))

        # 4. Verify GFF3 is parseable (if it exists)
        checks.append(self._check_gff_validity())

        # 5. Check for agent errors
        checks.append(self._check_agent_errors())

        passed = all(c["passed"] for c in checks)

        result = {
            "checks": checks,
            "all_passed": passed,
        }

        self.context.store_result(self.name, result)

        if passed:
            self.log.info("  Pipeline integrity [green]VERIFIED[/]")
        else:
            failed = [c["name"] for c in checks if not c["passed"]]
            self.log.warning(f"  Pipeline integrity issues: {', '.join(failed)}")
            self.send(
                "scrutinizer",
                MessageType.ALERT,
                {"issue": f"Pipeline integrity check failed: {', '.join(failed)}",
                 "severity": "warning"},
            )

        self.set_status(AgentStatus.DONE)
        return result

    def _check_locus_tags(self, genome) -> dict:
        """Ensure all locus tags are unique."""
        tags = [f.locus_tag for f in genome.all_features if f.locus_tag]
        unique = set(tags)
        duplicates = len(tags) - len(unique)

        if duplicates > 0:
            return {
                "name": "locus_tag_uniqueness",
                "passed": False,
                "message": f"{duplicates} duplicate locus tags found",
            }
        return {
            "name": "locus_tag_uniqueness",
            "passed": True,
            "message": f"All {len(tags)} locus tags are unique",
        }

    def _check_coordinates(self, genome) -> dict:
        """Verify all feature coordinates are within contig bounds."""
        out_of_bounds = 0
        for contig in genome.contigs:
            for f in contig.features:
                if f.start < 1 or f.end > contig.length:
                    out_of_bounds += 1
                if f.start > f.end:
                    out_of_bounds += 1

        if out_of_bounds > 0:
            return {
                "name": "coordinate_validity",
                "passed": False,
                "message": f"{out_of_bounds} features with invalid coordinates",
            }
        return {
            "name": "coordinate_validity",
            "passed": True,
            "message": "All feature coordinates are valid",
        }

    def _check_output_files(self, synth_result: dict) -> dict:
        """Check that output files exist and are non-empty."""
        files = synth_result.get("files_written", [])
        missing = []
        empty = []

        for f in files:
            path = Path(f)
            if not path.exists():
                missing.append(path.name)
            elif path.stat().st_size == 0:
                empty.append(path.name)

        if missing or empty:
            msg_parts = []
            if missing:
                msg_parts.append(f"missing: {', '.join(missing)}")
            if empty:
                msg_parts.append(f"empty: {', '.join(empty)}")
            return {
                "name": "output_files",
                "passed": False,
                "message": "; ".join(msg_parts),
            }
        return {
            "name": "output_files",
            "passed": True,
            "message": f"All {len(files)} output files verified",
        }

    def _check_gff_validity(self) -> dict:
        """Basic GFF3 format validation."""
        config = self.context.config
        gff_path = config.output_dir / f"{self.context.genome.name}.gff3"

        if not gff_path.exists():
            return {"name": "gff3_validity", "passed": True,
                    "message": "No GFF3 file to check (format not requested)"}

        try:
            with open(gff_path) as fh:
                first_line = fh.readline().strip()
                if not first_line.startswith("##gff-version"):
                    return {"name": "gff3_validity", "passed": False,
                            "message": "GFF3 missing version header"}

                # Check a few lines for tab-delimited format
                line_count = 0
                for line in fh:
                    if line.startswith("#"):
                        continue
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) != 9:
                        return {"name": "gff3_validity", "passed": False,
                                "message": f"GFF3 line has {len(parts)} columns (expected 9)"}
                    line_count += 1
                    if line_count >= 10:  # spot-check first 10
                        break

            return {"name": "gff3_validity", "passed": True,
                    "message": "GFF3 format valid"}

        except Exception as e:
            return {"name": "gff3_validity", "passed": False,
                    "message": f"GFF3 parse error: {e}"}

    def _check_agent_errors(self) -> dict:
        """Check if any agents recorded errors."""
        errors = self.context.errors
        if errors:
            return {
                "name": "agent_errors",
                "passed": False,
                "message": f"Agents with errors: {', '.join(errors.keys())}",
            }
        return {
            "name": "agent_errors",
            "passed": True,
            "message": "No agent errors recorded",
        }
