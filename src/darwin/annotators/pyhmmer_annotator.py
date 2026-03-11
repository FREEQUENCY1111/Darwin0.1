"""pyhmmer — protein function annotation via HMM profiles.

Uses pyhmmer (Python bindings for HMMER3) to search predicted proteins
against HMM profile databases (Pfam, TIGRFAMs, etc.) for functional
annotation. Much faster than subprocess calls to hmmsearch.

This annotator:
  1. Extracts protein sequences from CDS features
  2. Searches against provided HMM databases
  3. Assigns product names and domain info to features
"""

from __future__ import annotations

from pathlib import Path

import pyhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
from pyhmmer.plan7 import HMMFile

from darwin.annotators.base import BaseAnnotator
from darwin.models import AnnotationConfig, FeatureType, Genome
from darwin.utils.logging import get_logger


class PyhmmerAnnotator(BaseAnnotator):
    name = "pyhmmer"

    # Default HMM database search order
    DEFAULT_DBS = ["Pfam-A.hmm", "TIGRFAMs.hmm"]

    def __init__(self, config: AnnotationConfig, db_paths: list[Path] | None = None) -> None:
        super().__init__(config)
        self.db_paths = db_paths or []

    def check_dependencies(self) -> bool:
        """pyhmmer is a Python library — always available if installed."""
        try:
            import pyhmmer  # noqa: F401
            return True
        except ImportError:
            return False

    def run(self, genome: Genome) -> Genome:
        self.log.info("[bold blue]Running pyhmmer protein annotation...[/]")

        if not self.db_paths:
            self.log.warning("No HMM databases provided — skipping functional annotation")
            self.log.info(
                "  Download Pfam: wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
            )
            return genome

        # Collect all CDS features and their translated sequences
        proteins = self._extract_proteins(genome)
        if not proteins:
            self.log.warning("No CDS features found — run Prodigal first")
            return genome

        self.log.info(f"  Searching {len(proteins)} proteins against {len(self.db_paths)} HMM database(s)")

        # Convert to pyhmmer digital sequences
        alphabet = Alphabet.amino()
        sequences = DigitalSequenceBlock(alphabet)
        for name, seq in proteins:
            text_seq = TextSequence(name=name.encode(), sequence=seq)
            sequences.append(text_seq.digitize(alphabet))

        # Search each database
        all_hits: dict[str, tuple[str, float, str]] = {}  # locus_tag -> (product, evalue, db)

        for db_path in self.db_paths:
            self.log.info(f"  Scanning [cyan]{db_path.name}[/]...")
            hits = self._search_db(sequences, db_path)
            # Keep best hit per protein (lowest e-value)
            for locus_tag, product, evalue in hits:
                if locus_tag not in all_hits or evalue < all_hits[locus_tag][1]:
                    all_hits[locus_tag] = (product, evalue, db_path.name)

        # Apply annotations back to genome
        annotated = 0
        for contig in genome.contigs:
            for feature in contig.features:
                if feature.feature_type != FeatureType.CDS:
                    continue
                tag = feature.locus_tag
                if tag in all_hits:
                    product, evalue, db_name = all_hits[tag]
                    feature.attributes["product"] = product
                    feature.attributes["evalue"] = f"{evalue:.2e}"
                    feature.attributes["db_source"] = db_name
                    annotated += 1

        hypo = sum(
            1
            for c in genome.contigs
            for f in c.features
            if f.feature_type == FeatureType.CDS and f.product == "hypothetical protein"
        )
        self.log.info(
            f"  Annotated [green]{annotated}[/] proteins, "
            f"[yellow]{hypo}[/] remain hypothetical"
        )
        return genome

    def _extract_proteins(self, genome: Genome) -> list[tuple[str, str]]:
        """Extract (locus_tag, amino_acid_sequence) for all CDS features.

        Uses Biopython to translate nucleotide sequences from the genome.
        """
        from Bio.Seq import Seq

        proteins: list[tuple[str, str]] = []
        table = self.config.translation_table

        for contig in genome.contigs:
            for feature in contig.features:
                if feature.feature_type != FeatureType.CDS:
                    continue

                # Extract nucleotide sequence
                start = feature.start - 1  # convert to 0-based
                end = feature.end
                nuc_seq = contig.sequence[start:end]

                if feature.strand.value == "-":
                    nuc_seq = str(Seq(nuc_seq).reverse_complement())

                # Translate
                try:
                    protein = str(Seq(nuc_seq).translate(table=table, to_stop=True))
                    if len(protein) >= 20:  # skip very short ORFs
                        proteins.append((feature.locus_tag, protein))
                except Exception:
                    continue

        return proteins

    def _search_db(
        self,
        sequences: DigitalSequenceBlock,
        db_path: Path,
    ) -> list[tuple[str, str, float]]:
        """Search sequences against a single HMM database.

        Returns list of (locus_tag, product_name, evalue).
        """
        hits: list[tuple[str, str, float]] = []

        with HMMFile(str(db_path)) as hmm_file:
            for top_hits in pyhmmer.hmmsearch(
                hmm_file,
                sequences,
                cpus=self.config.cpus,
                E=self.config.evalue,
            ):
                for hit in top_hits:
                    if hit.included:
                        locus_tag = hit.name.decode()
                        # Use HMM description as product name
                        product = (
                            top_hits.query_name.decode()
                            if top_hits.query_name
                            else "hypothetical protein"
                        )
                        if top_hits.query_accession:
                            accession = top_hits.query_accession.decode()
                        else:
                            accession = ""

                        # Clean up product name
                        desc = ""
                        if hasattr(top_hits, "description") and top_hits.description:
                            desc = top_hits.description.decode()

                        final_product = desc if desc else product
                        hits.append((locus_tag, final_product, hit.evalue))

        return hits
