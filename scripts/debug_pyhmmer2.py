#!/usr/bin/env python3
"""Find correct pyhmmer TopHits attributes."""
import pyhmmer

print(f"pyhmmer version: {pyhmmer.__version__}")

alphabet = pyhmmer.easel.Alphabet.amino()
seq = pyhmmer.easel.TextSequence(name=b"test", sequence="MKKLLVLSLF").digitize(alphabet)

with pyhmmer.plan7.HMMFile("databases/Pfam-A.hmm") as hmm_file:
    for i, top_hits in enumerate(pyhmmer.hmmsearch(hmm_file, [seq], E=10)):
        if i == 0:
            print(f"\nTopHits type: {type(top_hits)}")
            print(f"TopHits dir: {[a for a in dir(top_hits) if not a.startswith('_')]}")
            if hasattr(top_hits, 'query_name'):
                print(f"query_name: {top_hits.query_name}")
            if hasattr(top_hits, 'query'):
                print(f"query: {top_hits.query}")
            if hasattr(top_hits, 'qname'):
                print(f"qname: {top_hits.qname}")
            break
