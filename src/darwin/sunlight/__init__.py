"""
Sunlight — The only external energy source.

Sunlight is how the user interacts with the jar.
It comes in two forms:
  - CLI: direct terminal access (the sun)
  - API: HTTP endpoints (artificial light)

Both do the same thing: send energy into the jar.
"""

from darwin.sunlight.cli import cli
from darwin.sunlight.api import create_app

__all__ = ["cli", "create_app"]
