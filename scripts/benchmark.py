#!/usr/bin/env python3
"""
Darwin Benchmark — Compare Darwin annotations against NCBI RefSeq gold standard.

Parses both GFF3 files, extracts CDS/tRNA/rRNA features, and computes:
  - Sensitivity (recall): what fraction of NCBI genes did Darwin find?
  - Precision: what fraction of Darwin's calls match a real NCBI gene?
  - Exact matches: identical start + stop + strand
  - Partial overlaps: reciprocal overlap >= 50%
  - Novel calls: Darwin features with no NCBI match
  - Missed genes: NCBI features Darwin didn't find

Usage:
    python scripts/benchmark.py \\
        --ncbi ecoli_benchmark/ncbi_dataset/data/GCF_000005845.2/genomic.gff \\
        --darwin darwin_vs_ncbi/GCF_000005845.2_ASM584v2_genomic.gff3 \\
        --output benchmark_report.json
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class Feature:
    """A genomic feature with coordinates."""

    seqid: str
    source: str
    ftype: str
    start: int
    end: int
    strand: str
    attributes: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def overlap(self, other: Feature) -> int:
        """Calculate base-pair overlap between two features."""
        if self.seqid != other.seqid:
            return 0
        start = max(self.start, other.start)
        end = min(self.end, other.end)
        return max(0, end - start + 1)

    def reciprocal_overlap(self, other: Feature) -> float:
        """Fraction of the smaller feature covered by overlap."""
        ovl = self.overlap(other)
        if ovl == 0:
            return 0.0
        min_len = min(self.length, other.length)
        return ovl / min_len

    def is_exact_match(self, other: Feature) -> bool:
        """Exact coordinate + strand match."""
        return (
            self.seqid == other.seqid
            and self.start == other.start
            and self.end == other.end
            and self.strand == other.strand
        )

    def get_attr(self, key: str) -> str:
        """Extract a value from GFF3 attributes."""
        for part in self.attributes.split(";"):
            if part.startswith(f"{key}="):
                return part[len(key) + 1 :]
        return ""


@dataclass
class BenchmarkResult:
    """Results of comparing Darwin vs NCBI annotations."""

    # Counts
    ncbi_cds: int = 0
    darwin_cds: int = 0
    ncbi_trna: int = 0
    darwin_trna: int = 0
    ncbi_rrna: int = 0
    darwin_rrna: int = 0

    # CDS comparison
    cds_exact_matches: int = 0
    cds_partial_matches: int = 0  # >= 50% reciprocal overlap
    cds_missed: int = 0  # NCBI genes Darwin didn't find
    cds_novel: int = 0  # Darwin calls with no NCBI match

    # tRNA comparison
    trna_exact_matches: int = 0
    trna_partial_matches: int = 0
    trna_missed: int = 0
    trna_novel: int = 0

    # rRNA comparison
    rrna_exact_matches: int = 0
    rrna_partial_matches: int = 0
    rrna_missed: int = 0
    rrna_novel: int = 0

    # Metrics
    cds_sensitivity: float = 0.0  # (exact + partial) / ncbi_cds
    cds_precision: float = 0.0  # (exact + partial) / darwin_cds
    cds_exact_rate: float = 0.0  # exact / ncbi_cds

    # Start codon accuracy
    start_codon_exact: int = 0
    start_codon_within_15bp: int = 0
    start_codon_offset_distribution: dict = field(default_factory=dict)

    # Missed gene details
    missed_gene_names: list = field(default_factory=list)
    novel_darwin_ids: list = field(default_factory=list)


def parse_gff3(path: Path, source_filter: str | None = None) -> list[Feature]:
    """Parse a GFF3 file and return features."""
    features = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if source_filter and parts[1] != source_filter:
                continue
            feat = Feature(
                seqid=parts[0],
                source=parts[1],
                ftype=parts[2],
                start=int(parts[3]),
                end=int(parts[4]),
                strand=parts[6],
                attributes=parts[8],
            )
            features.append(feat)
    return features


def extract_cds(features: list[Feature]) -> list[Feature]:
    """Extract CDS features."""
    return [f for f in features if f.ftype == "CDS"]


def extract_trna(features: list[Feature]) -> list[Feature]:
    """Extract tRNA features."""
    return [f for f in features if f.ftype == "tRNA"]


def extract_rrna(features: list[Feature]) -> list[Feature]:
    """Extract rRNA features."""
    return [f for f in features if f.ftype == "rRNA"]


def find_best_match(
    query: Feature, targets: list[Feature], min_overlap: float = 0.5
) -> tuple[Feature | None, str]:
    """
    Find the best matching target for a query feature.
    Returns (best_match, match_type) where match_type is 'exact', 'partial', or 'none'.
    """
    best = None
    best_ovl = 0.0

    for target in targets:
        if query.is_exact_match(target):
            return target, "exact"

        ovl = query.reciprocal_overlap(target)
        if ovl > best_ovl:
            best_ovl = ovl
            best = target

    if best and best_ovl >= min_overlap:
        return best, "partial"

    return None, "none"


def compare_feature_sets(
    ncbi_features: list[Feature],
    darwin_features: list[Feature],
    min_overlap: float = 0.5,
) -> dict:
    """Compare two sets of features and return match statistics."""
    exact = 0
    partial = 0
    missed = 0
    novel = 0
    missed_details = []
    novel_details = []

    # Start codon analysis (for CDS)
    start_exact = 0
    start_within_15 = 0
    start_offsets = {}

    # Track which Darwin features got matched
    matched_darwin = set()

    # For each NCBI feature, find best Darwin match
    for ncbi_feat in ncbi_features:
        match, mtype = find_best_match(ncbi_feat, darwin_features, min_overlap)

        if mtype == "exact":
            exact += 1
            start_exact += 1
            if match:
                matched_darwin.add(id(match))
        elif mtype == "partial":
            partial += 1
            if match:
                matched_darwin.add(id(match))
                # Analyze start codon offset
                if ncbi_feat.strand == match.strand:
                    if ncbi_feat.strand == "+":
                        offset = match.start - ncbi_feat.start
                    else:
                        offset = ncbi_feat.end - match.end
                    offset_key = str(offset)
                    start_offsets[offset_key] = start_offsets.get(offset_key, 0) + 1
                    if abs(offset) <= 15:
                        start_within_15 += 1
        else:
            missed += 1
            gene_name = ncbi_feat.get_attr("gene") or ncbi_feat.get_attr("Name") or ncbi_feat.get_attr("ID")
            missed_details.append(
                {
                    "name": gene_name,
                    "start": ncbi_feat.start,
                    "end": ncbi_feat.end,
                    "strand": ncbi_feat.strand,
                    "length": ncbi_feat.length,
                    "product": ncbi_feat.get_attr("product"),
                }
            )

    # Find novel Darwin calls (no NCBI match)
    for darwin_feat in darwin_features:
        if id(darwin_feat) not in matched_darwin:
            # Double-check: is there ANY overlap with NCBI?
            _, mtype = find_best_match(darwin_feat, ncbi_features, min_overlap)
            if mtype == "none":
                novel += 1
                novel_details.append(
                    {
                        "id": darwin_feat.get_attr("ID"),
                        "start": darwin_feat.start,
                        "end": darwin_feat.end,
                        "strand": darwin_feat.strand,
                        "length": darwin_feat.length,
                    }
                )

    return {
        "exact": exact,
        "partial": partial,
        "missed": missed,
        "novel": novel,
        "missed_details": missed_details[:20],  # Top 20
        "novel_details": novel_details[:20],
        "start_exact": start_exact,
        "start_within_15": start_within_15,
        "start_offsets": dict(sorted(start_offsets.items(), key=lambda x: abs(int(x[0])))[:20]),
    }


def run_benchmark(ncbi_path: Path, darwin_path: Path) -> BenchmarkResult:
    """Run the full benchmark comparison."""
    print("Parsing NCBI GFF3...")
    ncbi_all = parse_gff3(ncbi_path)  # No source filter — NCBI GFFs use varied sources
    ncbi_cds = extract_cds(ncbi_all)
    ncbi_trna = extract_trna(ncbi_all)
    ncbi_rrna = extract_rrna(ncbi_all)

    print("Parsing Darwin GFF3...")
    darwin_all = parse_gff3(darwin_path, source_filter="Darwin")
    darwin_cds = extract_cds(darwin_all)
    darwin_trna = extract_trna(darwin_all)
    darwin_rrna = extract_rrna(darwin_all)

    result = BenchmarkResult()
    result.ncbi_cds = len(ncbi_cds)
    result.darwin_cds = len(darwin_cds)
    result.ncbi_trna = len(ncbi_trna)
    result.darwin_trna = len(darwin_trna)
    result.ncbi_rrna = len(ncbi_rrna)
    result.darwin_rrna = len(darwin_rrna)

    print(f"\nNCBI:   {result.ncbi_cds} CDS, {result.ncbi_trna} tRNA, {result.ncbi_rrna} rRNA")
    print(f"Darwin: {result.darwin_cds} CDS, {result.darwin_trna} tRNA, {result.darwin_rrna} rRNA")

    # Compare CDS
    print("\nComparing CDS features...")
    cds_result = compare_feature_sets(ncbi_cds, darwin_cds)
    result.cds_exact_matches = cds_result["exact"]
    result.cds_partial_matches = cds_result["partial"]
    result.cds_missed = cds_result["missed"]
    result.cds_novel = cds_result["novel"]
    result.start_codon_exact = cds_result["start_exact"]
    result.start_codon_within_15bp = cds_result["start_within_15"]
    result.start_codon_offset_distribution = cds_result["start_offsets"]
    result.missed_gene_names = [m["name"] for m in cds_result["missed_details"]]
    result.novel_darwin_ids = [n["id"] for n in cds_result["novel_details"]]

    # Metrics
    total_found = result.cds_exact_matches + result.cds_partial_matches
    if result.ncbi_cds > 0:
        result.cds_sensitivity = round(total_found / result.ncbi_cds * 100, 2)
        result.cds_exact_rate = round(result.cds_exact_matches / result.ncbi_cds * 100, 2)
    if result.darwin_cds > 0:
        result.cds_precision = round(total_found / result.darwin_cds * 100, 2)

    # Compare tRNA
    print("Comparing tRNA features...")
    trna_result = compare_feature_sets(ncbi_trna, darwin_trna)
    result.trna_exact_matches = trna_result["exact"]
    result.trna_partial_matches = trna_result["partial"]
    result.trna_missed = trna_result["missed"]
    result.trna_novel = trna_result["novel"]

    # Compare rRNA
    print("Comparing rRNA features...")
    rrna_result = compare_feature_sets(ncbi_rrna, darwin_rrna)
    result.rrna_exact_matches = rrna_result["exact"]
    result.rrna_partial_matches = rrna_result["partial"]
    result.rrna_missed = rrna_result["missed"]
    result.rrna_novel = rrna_result["novel"]

    return result


def print_report(result: BenchmarkResult, ncbi_path: Path | None = None) -> None:
    """Print a human-readable benchmark report."""
    # Derive genome accession from NCBI GFF path (e.g. .../GCF_000005845.2/genomic.gff)
    label = "Unknown Genome"
    if ncbi_path:
        for part in ncbi_path.parts:
            if part.startswith("GCF_") or part.startswith("GCA_"):
                label = part
                break
        else:
            label = ncbi_path.stem

    print("\n" + "=" * 60)
    print("  DARWIN vs NCBI RefSeq — Benchmark Report")
    print(f"  {label}")
    print("=" * 60)

    print("\n--- Feature Counts ---")
    print(f"{'':>20} {'NCBI':>8} {'Darwin':>8}")
    print(f"{'CDS':>20} {result.ncbi_cds:>8} {result.darwin_cds:>8}")
    print(f"{'tRNA':>20} {result.ncbi_trna:>8} {result.darwin_trna:>8}")
    print(f"{'rRNA':>20} {result.ncbi_rrna:>8} {result.darwin_rrna:>8}")

    print("\n--- CDS Accuracy ---")
    print(f"  Exact matches (start+stop+strand):  {result.cds_exact_matches}")
    print(f"  Partial matches (>=50% overlap):     {result.cds_partial_matches}")
    total = result.cds_exact_matches + result.cds_partial_matches
    print(f"  Total matched:                       {total}")
    print(f"  Missed NCBI genes:                   {result.cds_missed}")
    print(f"  Novel Darwin calls:                  {result.cds_novel}")
    print(f"")
    print(f"  Sensitivity (recall):  {result.cds_sensitivity}%")
    print(f"  Precision:             {result.cds_precision}%")
    print(f"  Exact match rate:      {result.cds_exact_rate}%")

    print("\n--- Start Codon Analysis ---")
    print(f"  Exact start position:     {result.start_codon_exact}")
    print(f"  Within 15bp of NCBI:      {result.start_codon_within_15bp}")
    if result.start_codon_offset_distribution:
        print(f"  Top start offsets (bp):   ", end="")
        top = list(result.start_codon_offset_distribution.items())[:10]
        print(", ".join(f"{k}bp:{v}" for k, v in top))

    print("\n--- tRNA Accuracy ---")
    print(f"  Exact matches:  {result.trna_exact_matches}/{result.ncbi_trna}")
    print(f"  Partial:        {result.trna_partial_matches}")
    print(f"  Missed:         {result.trna_missed}")
    print(f"  Novel:          {result.trna_novel}")

    print("\n--- rRNA Accuracy ---")
    print(f"  Exact matches:  {result.rrna_exact_matches}/{result.ncbi_rrna}")
    print(f"  Partial:        {result.rrna_partial_matches}")
    print(f"  Missed:         {result.rrna_missed}")
    print(f"  Novel:          {result.rrna_novel}")

    if result.missed_gene_names:
        print(f"\n--- Sample Missed Genes (first 20) ---")
        for name in result.missed_gene_names[:20]:
            print(f"  • {name}")

    if result.novel_darwin_ids:
        print(f"\n--- Sample Novel Darwin Calls (first 20) ---")
        for did in result.novel_darwin_ids[:20]:
            print(f"  • {did}")

    # Overall grade
    print("\n" + "=" * 60)
    if result.cds_sensitivity >= 95:
        grade = "A — Excellent"
    elif result.cds_sensitivity >= 90:
        grade = "B — Good"
    elif result.cds_sensitivity >= 80:
        grade = "C — Acceptable"
    elif result.cds_sensitivity >= 70:
        grade = "D — Needs Work"
    else:
        grade = "F — Poor"
    print(f"  Overall CDS Grade: {grade}")
    print("=" * 60)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark Darwin vs NCBI RefSeq")
    parser.add_argument("--ncbi", type=Path, required=True, help="NCBI reference GFF3")
    parser.add_argument("--darwin", type=Path, required=True, help="Darwin output GFF3")
    parser.add_argument("--output", type=Path, default=None, help="Output JSON report")
    args = parser.parse_args()

    if not args.ncbi.exists():
        print(f"Error: NCBI GFF3 not found: {args.ncbi}", file=sys.stderr)
        sys.exit(1)
    if not args.darwin.exists():
        print(f"Error: Darwin GFF3 not found: {args.darwin}", file=sys.stderr)
        sys.exit(1)

    result = run_benchmark(args.ncbi, args.darwin)
    print_report(result, ncbi_path=args.ncbi)

    if args.output:
        with open(args.output, "w") as f:
            json.dump(asdict(result), f, indent=2)
        print(f"\nJSON report saved to {args.output}")


if __name__ == "__main__":
    main()
