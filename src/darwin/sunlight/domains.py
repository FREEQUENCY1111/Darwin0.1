"""
Domain characterization engine — standalone protein domain analysis.

Takes a protein FASTA and searches it against Pfam HMM databases,
returning all domain hits per protein with positions, scores, and
an ASCII architecture diagram.

Unlike the genome annotation pipeline (which keeps only the best hit
per protein), this returns ALL domain hits so you can see the full
domain architecture of each protein.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger("darwin.domains")


@dataclass
class DomainHit:
    """A single domain hit on a protein."""

    protein_id: str
    protein_length: int
    domain_name: str
    accession: str
    description: str
    env_from: int  # envelope start (1-based)
    env_to: int    # envelope end (1-based)
    evalue: float
    score: float
    db_name: str

    @property
    def domain_length(self) -> int:
        return self.env_to - self.env_from + 1

    @property
    def coverage(self) -> float:
        """Fraction of the protein covered by this domain."""
        if self.protein_length == 0:
            return 0.0
        return self.domain_length / self.protein_length


@dataclass
class ProteinDomains:
    """All domains found on a single protein."""

    protein_id: str
    protein_length: int
    hits: list[DomainHit] = field(default_factory=list)

    @property
    def num_domains(self) -> int:
        return len(self.hits)

    @property
    def total_coverage(self) -> float:
        """Fraction of protein covered by any domain (merging overlaps)."""
        if not self.hits or self.protein_length == 0:
            return 0.0
        # Merge overlapping intervals
        intervals = sorted((h.env_from, h.env_to) for h in self.hits)
        merged: list[tuple[int, int]] = [intervals[0]]
        for start, end in intervals[1:]:
            if start <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        covered = sum(e - s + 1 for s, e in merged)
        return covered / self.protein_length


def read_protein_fasta(fasta_path: Path) -> dict[str, str]:
    """Read a protein FASTA file and return {id: sequence} dict."""
    proteins: dict[str, str] = {}
    current_id = ""
    current_seq: list[str] = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id:
                    proteins[current_id] = "".join(current_seq)
                # Take first word after > as ID
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            proteins[current_id] = "".join(current_seq)

    return proteins


def search_domains(
    proteins: dict[str, str],
    hmm_paths: list[Path],
    evalue: float = 1e-5,
    cpus: int = 1,
) -> list[ProteinDomains]:
    """
    Search proteins against HMM database(s) and return all domain hits.

    Unlike the genome pipeline's PyhmmerPlant (which keeps only the best
    hit per protein), this collects ALL domain hits with envelope coordinates
    so we can map the full domain architecture.
    """
    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = []
    length_map: dict[str, int] = {}
    tag_map: dict[str | bytes, str] = {}

    for tag, seq in proteins.items():
        try:
            ds = pyhmmer.easel.TextSequence(
                name=tag.encode(),
                sequence=seq,
            ).digitize(alphabet)
            sequences.append(ds)
            length_map[tag] = len(seq)
            tag_map[tag] = tag
            tag_map[tag.encode()] = tag
        except Exception:
            logger.warning(f"Skipping invalid sequence: {tag}")
            continue

    if not sequences:
        return []

    # Collect all domain hits per protein
    hits_by_protein: dict[str, list[DomainHit]] = {}

    for hmm_path in hmm_paths:
        db_name = hmm_path.stem
        logger.info(f"Searching {len(sequences)} proteins against {db_name}...")

        with pyhmmer.plan7.HMMFile(str(hmm_path)) as hmm_file:
            for top_hits in pyhmmer.hmmsearch(
                hmm_file, sequences, E=evalue, cpus=cpus
            ):
                # top_hits.query is the HMM profile
                qname = top_hits.query.name
                hmm_name = qname.decode() if isinstance(qname, bytes) else qname

                qacc = getattr(top_hits.query, "accession", None)
                if qacc:
                    acc = qacc.decode() if isinstance(qacc, bytes) else qacc
                else:
                    acc = hmm_name

                qdesc = getattr(top_hits.query, "description", None)
                if qdesc:
                    desc = qdesc.decode() if isinstance(qdesc, bytes) else qdesc
                else:
                    desc = ""

                for hit in top_hits:
                    if not hit.included:
                        continue

                    tag = tag_map.get(hit.name, "")
                    if not tag:
                        continue

                    # Iterate over domain annotations for this hit
                    for domain in hit.domains:
                        if not domain.included:
                            continue

                        protein_len = length_map.get(tag, 0)
                        domain_hit = DomainHit(
                            protein_id=tag,
                            protein_length=protein_len,
                            domain_name=hmm_name,
                            accession=acc,
                            description=desc,
                            env_from=domain.env_from,
                            env_to=domain.env_to,
                            evalue=domain.i_evalue,
                            score=domain.score,
                            db_name=db_name,
                        )

                        if tag not in hits_by_protein:
                            hits_by_protein[tag] = []
                        hits_by_protein[tag].append(domain_hit)

    # Build ProteinDomains results (sorted by protein ID)
    results: list[ProteinDomains] = []
    for tag in sorted(proteins.keys()):
        hits = hits_by_protein.get(tag, [])
        # Sort domains by position on the protein
        hits.sort(key=lambda h: (h.env_from, h.env_to))
        results.append(
            ProteinDomains(
                protein_id=tag,
                protein_length=length_map.get(tag, len(proteins[tag])),
                hits=hits,
            )
        )

    return results


def format_tsv(results: list[ProteinDomains]) -> str:
    """Format domain results as a TSV table."""
    header = "\t".join([
        "protein_id",
        "protein_length",
        "domain_name",
        "accession",
        "description",
        "env_from",
        "env_to",
        "domain_length",
        "evalue",
        "score",
        "coverage",
    ])
    lines = [header]

    for prot in results:
        if not prot.hits:
            # Still list the protein with no domains
            lines.append("\t".join([
                prot.protein_id,
                str(prot.protein_length),
                "-",
                "-",
                "no domains found",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
            ]))
        else:
            for hit in prot.hits:
                lines.append("\t".join([
                    hit.protein_id,
                    str(hit.protein_length),
                    hit.domain_name,
                    hit.accession,
                    hit.description,
                    str(hit.env_from),
                    str(hit.env_to),
                    str(hit.domain_length),
                    f"{hit.evalue:.2e}",
                    f"{hit.score:.1f}",
                    f"{hit.coverage:.1%}",
                ]))

    return "\n".join(lines) + "\n"


def format_domain_map(prot: ProteinDomains, width: int = 70) -> str:
    """
    Generate an ASCII domain architecture map for a single protein.

    Example output:
        RecA (354 aa, 2 domains, 85.3% covered)
        |=====[AAA_ATPase]=========------[RecA_C]==|
        1                                         354
    """
    if prot.protein_length == 0:
        return f"  {prot.protein_id} (0 aa)\n"

    label = prot.protein_id
    n_domains = prot.num_domains
    cov = prot.total_coverage

    header = f"  {label} ({prot.protein_length} aa"
    if n_domains > 0:
        header += f", {n_domains} domain{'s' if n_domains != 1 else ''}, {cov:.1%} covered"
    else:
        header += ", no domains found"
    header += ")"

    if n_domains == 0:
        bar = "|" + "-" * (width - 2) + "|"
        positions = f"  1{' ' * (width - len(str(prot.protein_length)) - 1)}{prot.protein_length}"
        return f"{header}\n  {bar}\n{positions}\n"

    # Build the bar character by character
    scale = (width - 2) / prot.protein_length  # -2 for the | delimiters
    bar = ["-"] * (width - 2)

    # Place each domain on the bar
    domain_labels: list[tuple[int, int, str]] = []
    for hit in prot.hits:
        col_start = int((hit.env_from - 1) * scale)
        col_end = int((hit.env_to - 1) * scale)
        col_start = max(0, min(col_start, width - 3))
        col_end = max(col_start, min(col_end, width - 3))

        # Fill with = for domain region
        for i in range(col_start, col_end + 1):
            bar[i] = "="

        # Try to place a label
        short_name = hit.domain_name
        if len(short_name) > 15:
            short_name = short_name[:14] + "~"
        label_text = f"[{short_name}]"
        label_len = len(label_text)

        # Center the label in the domain region
        region_len = col_end - col_start + 1
        if region_len >= label_len:
            label_start = col_start + (region_len - label_len) // 2
            domain_labels.append((label_start, label_len, label_text))

    # Apply labels to the bar
    for start, length, text in domain_labels:
        for i, ch in enumerate(text):
            if start + i < len(bar):
                bar[start + i] = ch

    bar_str = "|" + "".join(bar) + "|"
    positions = f"  1{' ' * (width - len(str(prot.protein_length)) - 1)}{prot.protein_length}"

    return f"{header}\n  {bar_str}\n{positions}\n"


def format_all_maps(results: list[ProteinDomains], width: int = 70) -> str:
    """Generate ASCII domain maps for all proteins."""
    lines: list[str] = []
    lines.append("=" * (width + 4))
    lines.append("  DOMAIN ARCHITECTURE MAP")
    lines.append("=" * (width + 4))
    lines.append("")

    for prot in results:
        lines.append(format_domain_map(prot, width=width))

    # Legend
    lines.append("-" * (width + 4))
    lines.append("  Legend:  = domain region   - uncharacterized   | protein boundary")
    lines.append("")

    return "\n".join(lines)
