"""Base class for all annotators."""

from __future__ import annotations

import abc
from pathlib import Path

from darwin.models import AnnotationConfig, Genome
from darwin.utils.logging import get_logger


class BaseAnnotator(abc.ABC):
    """Abstract annotator that all tools must implement."""

    name: str = "base"

    def __init__(self, config: AnnotationConfig) -> None:
        self.config = config
        self.log = get_logger(f"darwin.{self.name}")
        self._work_dir = config.output_dir / ".darwin_tmp" / self.name
        self._work_dir.mkdir(parents=True, exist_ok=True)

    @abc.abstractmethod
    def run(self, genome: Genome) -> Genome:
        """Run annotation and return the genome with new features added.

        Implementations should:
          1. Write temp input files to self._work_dir
          2. Run the tool
          3. Parse output and add Features to genome.contigs[*].features
          4. Return the modified genome
        """
        ...

    @abc.abstractmethod
    def check_dependencies(self) -> bool:
        """Verify that all required tools / databases are available."""
        ...

    @property
    def work_dir(self) -> Path:
        return self._work_dir
