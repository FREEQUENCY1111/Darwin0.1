"""Synthesizer Agent — combines all annotations into final output files.

Responsibilities:
  - Collect the annotated genome from context
  - Generate all requested output formats (GFF3, GenBank, JSON, FASTA)
  - Create a summary report
  - Compile the full audit trail
"""

from __future__ import annotations

from typing import Any

from Bio.Seq import Seq

from darwin.agents.base import Agent, AgentStatus, MessageType
from darwin.output.gff import write_gff3
from darwin.output.genbank import write_genbank
from darwin.output.json_out import write_json
from darwin.utils.fasta import write_fasta, write_proteins


class Synthesizer(Agent):
    name = "synthesizer"
    role = "Combines all annotations into final output files and reports"

    def validate_input(self, **kwargs: Any) -> bool:
        return self.context.genome is not None

    def execute(self, **kwargs: Any) -> dict[str, Any]:
        self.set_status(AgentStatus.WORKING)
        self.log.info("[bold cyan]Synthesizer[/] — generating output files")

        genome = self.context.get_genome()
        config = self.context.config
        out = config.output_dir
        out.mkdir(parents=True, exist_ok=True)
        name = genome.name

        formats = config.formats
        written: list[str] = []

        if "gff3" in formats:
            path = write_gff3(genome, out / f"{name}.gff3")
            written.append(str(path))

        if "gbk" in formats:
            path = write_genbank(genome, out / f"{name}.gbk")
            written.append(str(path))

        if "json" in formats:
            path = write_json(genome, out / f"{name}.json")
            written.append(str(path))

        if "fna" in formats:
            path = write_fasta(genome.contigs, out / f"{name}.fna")
            written.append(str(path))

        if "faa" in formats:
            proteins = self._extract_proteins(genome)
            if proteins:
                path = write_proteins(proteins, out / f"{name}.faa")
                written.append(str(path))

        # Write audit trail
        audit = self.context.get_full_audit()
        import json
        audit_path = out / f"{name}_audit.json"
        with open(audit_path, "w") as fh:
            json.dump(audit, fh, indent=2, default=str)
        written.append(str(audit_path))

        result = {
            "files_written": written,
            "formats": formats,
            "genome_summary": genome.summary(),
        }

        self.context.store_result(self.name, result)

        self.broadcast(
            MessageType.STATUS,
            {"message": f"Synthesizer wrote {len(written)} output files"},
        )

        for w in written:
            self.log.info(f"  Wrote: [cyan]{w}[/]")

        self.set_status(AgentStatus.DONE)
        return result

    def _extract_proteins(self, genome) -> list[tuple[str, str, str]]:
        """Extract protein sequences from CDS features."""
        from darwin.models import FeatureType

        proteins = []
        table = self.context.config.translation_table

        for contig in genome.contigs:
            for f in contig.features:
                if f.feature_type != FeatureType.CDS:
                    continue
                nuc = contig.sequence[f.start - 1 : f.end]
                if f.strand.value == "-":
                    nuc = str(Seq(nuc).reverse_complement())
                try:
                    aa = str(Seq(nuc).translate(table=table, to_stop=True))
                    proteins.append((f.locus_tag, f.product, aa))
                except Exception:
                    continue
        return proteins
