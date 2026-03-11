"""Helpers for running external bioinformatics tools."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from darwin.utils.logging import get_logger

log = get_logger(__name__)


class ToolNotFoundError(Exception):
    """Raised when a required external tool is not on PATH."""


def which_tool(name: str) -> Path:
    """Find an external tool on PATH, or raise a clear error."""
    found = shutil.which(name)
    if found is None:
        raise ToolNotFoundError(
            f"Required tool '{name}' not found on PATH. "
            f"Install it with: sudo apt-get install {name}  (or conda install -c bioconda {name})"
        )
    return Path(found)


def run_external(
    cmd: list[str],
    *,
    check: bool = True,
    capture: bool = True,
    timeout: int = 600,
) -> subprocess.CompletedProcess:
    """Run an external command with logging and error handling.

    Args:
        cmd: command and arguments as a list.
        check: raise on non-zero exit code.
        capture: capture stdout/stderr.
        timeout: seconds before killing the process.

    Returns:
        CompletedProcess result.
    """
    cmd_str = " ".join(str(c) for c in cmd)
    log.debug(f"Running: {cmd_str}")
    try:
        result = subprocess.run(
            cmd,
            check=check,
            capture_output=capture,
            text=True,
            timeout=timeout,
        )
        return result
    except subprocess.TimeoutExpired:
        log.error(f"Command timed out after {timeout}s: {cmd_str}")
        raise
    except subprocess.CalledProcessError as exc:
        log.error(f"Command failed (exit {exc.returncode}): {cmd_str}")
        if exc.stderr:
            log.error(f"stderr: {exc.stderr[:500]}")
        raise
