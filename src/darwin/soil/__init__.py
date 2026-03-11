"""
Soil — The nutrient store.

Soil holds the reference databases, HMM profiles, and lookup tables
that feed the flora. Plants draw nutrients from soil to grow.

Without soil, flora can't produce meaningful annotations.
Soil doesn't act — it just stores and provides.
"""

from darwin.soil.nutrients import NutrientStore

__all__ = ["NutrientStore"]
