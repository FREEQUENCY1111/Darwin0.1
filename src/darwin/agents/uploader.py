"""Uploader Agent — exports and delivers results to the user.

Responsibilities:
  - Copy final outputs to the designated delivery location
  - Generate download links (for API mode)
  - Optionally push results to external services (NCBI, S3, etc.)
  - Create a summary notification
  - Clean up temporary working files
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

from darwin.agents.base import Agent, AgentStatus, MessageType


class Uploader(Agent):
    name = "uploader"
    role = "Exports and delivers final results to the user"

    def validate_input(self, **kwargs: Any) -> bool:
        synth_result = self.context.get_result("synthesizer")
        if not synth_result or not synth_result.get("files_written"):
            self.log.error("No output files — Synthesizer must run first")
            return False
        return True

    def execute(self, **kwargs: Any) -> dict[str, Any]:
        self.set_status(AgentStatus.WORKING)
        self.log.info("[bold magenta]Uploader[/] — delivering results")

        synth_result = self.context.get_result("synthesizer")
        files = synth_result.get("files_written", [])
        delivery_dir = kwargs.get("delivery_dir")

        delivered: list[str] = []

        if delivery_dir:
            # Copy files to delivery directory
            delivery_path = Path(delivery_dir)
            delivery_path.mkdir(parents=True, exist_ok=True)
            for f in files:
                src = Path(f)
                if src.exists():
                    dst = delivery_path / src.name
                    shutil.copy2(src, dst)
                    delivered.append(str(dst))
        else:
            # Files are already in output_dir
            delivered = files

        # Clean up temp working files
        tmp_dir = self.context.config.output_dir / ".darwin_tmp"
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir, ignore_errors=True)
            self.log.info("  Cleaned up temporary files")

        result = {
            "delivered_files": delivered,
            "delivery_location": str(delivery_dir or self.context.config.output_dir),
            "file_count": len(delivered),
        }

        self.context.store_result(self.name, result)

        self.broadcast(
            MessageType.STATUS,
            {"message": f"Results delivered: {len(delivered)} files"},
        )

        self.log.info(f"  Delivered [green]{len(delivered)}[/] files")
        self.set_status(AgentStatus.DONE)
        return result
