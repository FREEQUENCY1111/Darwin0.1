"""
NutrientStore — Database and reference data management.

Like minerals in soil, these are the raw materials that flora
draw from to produce annotations. The soil holds:
- HMM databases (Pfam, TIGRFAMs)
- Tool binaries and their locations
- Reference data for enrichment

Soil is passive — flora reaches into it when hungry.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger("darwin.soil")


@dataclass
class ToolInfo:
    """Information about an installed bioinformatics tool."""

    name: str
    binary: Path | None = None
    version: str = "unknown"
    available: bool = False

    def check(self) -> bool:
        """Check if this tool is planted in the soil."""
        path = shutil.which(self.name)
        if path:
            self.binary = Path(path)
            self.available = True
            logger.debug(f"🌱 {self.name} found at {path}")
        else:
            self.available = False
            logger.debug(f"🏜️ {self.name} not found in soil")
        return self.available


@dataclass
class HMMDatabase:
    """A single HMM database — concentrated nutrients."""

    name: str
    path: Path
    description: str = ""

    @property
    def available(self) -> bool:
        return self.path.exists()


class NutrientStore:
    """
    The soil layer — provides tools and databases to flora.

    Flora don't manage databases. They just reach into the soil
    and take what they need. If the soil is barren, the flora
    can't grow — but they handle that gracefully.
    """

    def __init__(self, hmm_databases: list[Path] | None = None) -> None:
        self._tools: dict[str, ToolInfo] = {}
        self._hmm_dbs: list[HMMDatabase] = []

        # Register core tools that flora might need
        for name in ["prodigal", "aragorn", "barrnap", "hmmsearch", "cmscan", "mob_recon", "isescan.py"]:
            self._tools[name] = ToolInfo(name=name)

        # Register explicitly provided HMM databases
        if hmm_databases:
            for db_path in hmm_databases:
                db_path = Path(db_path)
                self._hmm_dbs.append(
                    HMMDatabase(
                        name=db_path.stem,
                        path=db_path,
                    )
                )
        else:
            # Auto-discover cached databases from ~/.darwin/databases/
            self._discover_cached_databases()

    def survey(self) -> dict[str, bool]:
        """
        Survey the soil — what's available?

        Like testing soil quality before planting.
        Returns which tools and databases are present.
        """
        results = {}
        for name, tool in self._tools.items():
            tool.check()
            results[name] = tool.available

        for db in self._hmm_dbs:
            results[f"hmm:{db.name}"] = db.available

        available = sum(1 for v in results.values() if v)
        total = len(results)
        logger.info(f"🌱 Soil survey: {available}/{total} nutrients available")
        return results

    def get_tool(self, name: str) -> ToolInfo | None:
        """Reach into the soil for a specific tool."""
        return self._tools.get(name)

    def get_tool_path(self, name: str) -> Path | None:
        """Get the binary path for a tool, or None if not in soil."""
        tool = self._tools.get(name)
        if tool and tool.available:
            return tool.binary
        return None

    def get_hmm_databases(self) -> list[HMMDatabase]:
        """Get all available HMM databases."""
        return [db for db in self._hmm_dbs if db.available]

    @property
    def has_prodigal(self) -> bool:
        t = self._tools.get("prodigal")
        return bool(t and t.available)

    @property
    def has_aragorn(self) -> bool:
        t = self._tools.get("aragorn")
        return bool(t and t.available)

    @property
    def has_barrnap(self) -> bool:
        t = self._tools.get("barrnap")
        return bool(t and t.available)

    @property
    def has_mob_suite(self) -> bool:
        t = self._tools.get("mob_recon")
        return bool(t and t.available)

    @property
    def has_isescan(self) -> bool:
        t = self._tools.get("isescan.py")
        return bool(t and t.available)

    @property
    def has_hmm(self) -> bool:
        return len(self.get_hmm_databases()) > 0

    def _discover_cached_databases(self) -> None:
        """Auto-discover HMM databases from standard locations.

        Search order:
          1. ./databases/          (project-local)
          2. ~/.darwin/databases/  (user cache)
        """
        search_dirs = [
            Path.cwd() / "databases",
            Path.home() / ".darwin" / "databases",
        ]
        seen: set[str] = set()

        for search_dir in search_dirs:
            if not search_dir.exists():
                continue
            for hmm_file in sorted(search_dir.glob("*.hmm")):
                # De-duplicate by stem name (prefer first found)
                if hmm_file.stem in seen:
                    continue
                seen.add(hmm_file.stem)
                self._hmm_dbs.append(
                    HMMDatabase(
                        name=hmm_file.stem,
                        path=hmm_file,
                    )
                )
                logger.info(f"🌱 Auto-discovered DB: {hmm_file.name} ({search_dir})")

    @property
    def is_fertile(self) -> bool:
        """Is there enough in the soil for anything to grow?"""
        return self.has_prodigal  # at minimum, we need gene calling
