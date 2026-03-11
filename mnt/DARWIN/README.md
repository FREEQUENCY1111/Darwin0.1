# Darwin 🧬

Fast, accurate prokaryotic genome annotation.

Darwin takes a FASTA genome assembly and produces comprehensive annotations — protein-coding genes, tRNAs, rRNAs, and functional assignments — in seconds.

## Quick Start

```bash
# Install
pip install -e ".[dev,test]"

# Check external tools are available
darwin check

# Annotate a genome
darwin annotate genome.fasta -o results/

# Start the REST API
darwin serve --port 8000
```

## Annotation Pipeline

1. **Barrnap** — ribosomal RNA (5S, 16S, 23S) prediction
2. **Aragorn** — tRNA and tmRNA detection
3. **Prodigal** — ab initio protein-coding gene prediction
4. **pyhmmer** — HMM-based functional annotation (Pfam, TIGRFAMs)

## Output Formats

- GFF3 (`.gff3`)
- GenBank (`.gbk`)
- Protein FASTA (`.faa`)
- Nucleotide FASTA (`.fna`)
- JSON (`.json`)

## REST API

```bash
# Submit a genome
curl -X POST http://localhost:8000/annotate \
  -F "file=@genome.fasta"

# Check job status
curl http://localhost:8000/jobs/{job_id}
```

## External Dependencies

Install via conda: `conda install -c bioconda prodigal barrnap aragorn`

## License

MIT
