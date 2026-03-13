#!/usr/bin/env python3
"""Debug: verify pyhmmer 0.12 fix - str vs bytes."""

import pyhmmer

print(f"pyhmmer version: {pyhmmer.__version__}")

# Real E. coli dnaA protein (truncated but enough for Pfam hit)
test_proteins = {
    b"TEST_00002": (
        "MSLSLWQQCLARLQDELPATEFSMWIRPLQAELSDNTLALYAPNRFVLDWVRDKYLNNINGLLTSFCGADAPQL"
        "IRHFIDPESLHFSREIFASDSRRTPIAPATRRNFDRRPEFVNCIRTQKLEELINDAAATLTPLRHAALALLDSQ"
        "IGHGFFLTLRAQEELLSALKDKQFHRTRLHDAVSQTLIQHGIFVSALREKRSVLEAALETALDTISQLAALEH"
        "MGEFIDSLVEGALTKVYTEALELRTRGSITRLHLAQFDTFVELQKLTPQQVAALPQFMRDVLEATRHQISRHRA"
        "ALDLAARRGDHMEIKQFTEQMLQLTTRQHYQQYIADQFKLRESYQAFNIGNSAIDDNFHSLNLELYQAAQNLF"
    ),
}

alphabet = pyhmmer.easel.Alphabet.amino()
sequences = []
tag_map = {}

for name, seq in test_proteins.items():
    ds = pyhmmer.easel.TextSequence(name=name, sequence=seq).digitize(alphabet)
    sequences.append(ds)
    # Fixed: store both str and bytes keys
    tag_str = name.decode()
    tag_map[tag_str] = tag_str
    tag_map[name] = tag_str

print(f"tag_map keys: {list(tag_map.keys())}")

hit_count = 0
annotated = 0

try:
    with pyhmmer.plan7.HMMFile("databases/Pfam-A.hmm") as hmm_file:
        for hits in pyhmmer.hmmsearch(hmm_file, sequences, E=1e-10):
            for hit in hits:
                if hit.included:
                    hit_count += 1
                    tag = tag_map.get(hit.name, "")
                    qname = hits.query.name
                    product = qname.decode() if isinstance(qname, bytes) else qname

                    print(f"HIT: hit.name={hit.name!r} (type={type(hit.name).__name__})")
                    print(f"     tag_map lookup -> {tag!r}")
                    print(f"     query.name={qname!r} (type={type(qname).__name__})")
                    print(f"     product={product!r}")
                    print(f"     score={hit.score:.1f}")

                    if tag:
                        annotated += 1

    print(f"\nTotal hits: {hit_count}, Annotated: {annotated}")
    print("SUCCESS - fix works!" if annotated > 0 else "FAIL - still no annotations")

except Exception as e:
    print(f"\nEXCEPTION: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
