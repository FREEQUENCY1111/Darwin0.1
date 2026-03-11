"""Logging setup with Rich formatting."""

import logging
import sys
from rich.logging import RichHandler


def setup_logging(level: str = "INFO") -> None:
    """Configure logging with Rich output."""
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
    )


def get_logger(name: str) -> logging.Logger:
    """Get a named logger."""
    return logging.getLogger(name)
