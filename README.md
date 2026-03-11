# Darwin 🫙

**Prokaryotic Genome Annotator — A Self-Sustaining Ecosphere**

No orchestrators. No hierarchies. Just nature.

## The Metaphor

Picture a mason jar. Sterile. Empty. You go to the river, fill it with rocks, soil, water, plants, and organic matter. Close the lid. Set it in the sunlight. It becomes a self-sustaining ecosystem rooted in perfect harmony and balance.

That's Darwin.

## Architecture

| Layer | Metaphor | Code |
|-------|----------|------|
| 🫙 **Jar** | Sealed container | `darwin.jar` — runtime ecosphere |
| 🪨 **Rocks** | Immutable foundation | `darwin.rocks` — data models (Genome, Contig, Feature) |
| 💧 **Water** | Reactive stream | `darwin.water` — event bus connecting all life |
| 🌱 **Soil** | Nutrient store | `darwin.soil` — databases, tools, references |
| 🌿 **Flora** | Producers | `darwin.flora` — Prodigal, pyhmmer, Aragorn, Barrnap |
| 🦠 **Microbiome** | Decomposers | `darwin.microbiome` — QC, enrichment, output synthesis |
| ☀️ **Sunlight** | Energy input | `darwin.sunlight` — CLI and REST API |

## Quick Start

```bash
# Install
pip install -e ".[dev,test]"

# Check your soil (available tools)
darwin check

# Add sunlight (annotate a genome)
darwin annotate genome.fasta -o results/ -p MYORG

# Start the API server
darwin serve --port 8000
```

## How It Works

1. **Sunlight enters** — you submit a FASTA file
2. **Flora grow** — Prodigal, Aragorn, Barrnap all activate *simultaneously* (they all feed on `genome.loaded`)
3. **Epiphytes feed** — pyhmmer feeds on genes that Prodigal produced
4. **Microbiome decomposes** — Scrutinizer checks quality, Enricher adds context
5. **Synthesizer crystallizes** — output files are generated
6. **Equilibrium** — the jar is complete 🌊

No component knows about any other. They just react to nutrients in the water.

## API

```bash
# Submit a genome
curl -X POST http://localhost:8000/annotate \
  -F "file=@genome.fasta" \
  -F "locus_prefix=MYORG"

# Check job status
curl http://localhost:8000/jobs/{job_id}
```

## License

MIT
