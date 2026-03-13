#!/usr/bin/env python3
"""Quick diagnostic: are translations reaching pyhmmer? Are hits being found?"""

import sys
sys.path.insert(0, "src")

from pathlib import Path
from darwin.rocks.fasta import parse_fasta
from darwin.rocks.models import FeatureType
import asyncio
import tempfile

async def main():
    fasta = Path("ecoli_benchmark/ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna")
    genome = parse_fasta(fasta)
    print(f"Genome: {genome.name}, {genome.num_contigs} contigs")

    # Simulate Prodigal
    from darwin.flora.prodigal import ProdigalPlant
    from darwin.soil.nutrients import NutrientStore
    from darwin.water.stream import Stream, Nutrient, NutrientType

    stream = Stream()
    soil = NutrientStore(hmm_databases=[Path("databases/Pfam-A.hmm")])
    soil.survey()

    prodigal = ProdigalPlant(stream, soil)

    nutrient = Nutrient(
        type=NutrientType.GENOME_LOADED,
        data={"genome": genome, "config": {"locus_tag_prefix": "ECOLI"}},
        source="debug",
    )

    result = await prodigal.grow(nutrient)
    if not result:
        print("Prodigal failed!")
        return

    genome2 = result.data["genome"]
    cds = genome2.features_of_type(FeatureType.CDS)
    print(f"\nCDS features: {len(cds)}")

    with_translation = [f for f in cds if f.translation]
    print(f"With translation: {len(with_translation)}")
    without_translation = [f for f in cds if not f.translation]
    print(f"Without translation: {len(without_translation)}")

    if with_translation:
        f = with_translation[0]
        print(f"\nFirst translated CDS: {f.locus_tag}")
        print(f"  Translation (first 60 chars): {f.translation[:60]}")
    else:
        print("\n*** NO TRANSLATIONS FOUND ***")
        if cds:
            f = cds[0]
            print(f"  First CDS: {f.locus_tag}, translation='{f.translation}'")
        return

    # Now test pyhmmer matching
    print("\n--- Testing pyhmmer ---")
    try:
        import pyhmmer
        alphabet = pyhmmer.easel.Alphabet.amino()
        proteins = {f.locus_tag: f.translation for f in cds if f.translation}
        print(f"Proteins to search: {len(proteins)}")

        sequences = []
        tag_map = {}
        for tag, seq in proteins.items():
            try:
                ds = pyhmmer.easel.TextSequence(
                    name=tag.encode(),
                    sequence=seq,
                ).digitize(alphabet)
                sequences.append(ds)
                tag_map[tag.encode()] = tag
            except Exception as e:
                print(f"  Failed to digitize {tag}: {e}")

        print(f"Valid sequences: {len(sequences)}")
        print(f"First tag_map entry: {list(tag_map.items())[0]}")

        hmm_path = "databases/Pfam-A.hmm"
        print(f"\nSearching {hmm_path}...")
        hit_count = 0
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            for top_hits in pyhmmer.hmmsearch(hmm_file, sequences, E=1e-10):
                for hit in top_hits:
                    if hit.included:
                        hit_count += 1
                        if hit_count <= 5:
                            tag = tag_map.get(hit.name, "???")
                            print(f"  HIT: {top_hits.query_name.decode()} -> {tag} (score={hit.score:.1f})")
                        if hit_count == 6:
                            print("  ... (more hits)")

        print(f"\nTotal included hits: {hit_count}")
    except ImportError:
        print("pyhmmer not installed")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

asyncio.run(main())
