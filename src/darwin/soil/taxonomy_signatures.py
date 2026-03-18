"""
Taxonomy signatures — hierarchical probe-based classification
for major prokaryotic groups.

Uses diagnostic oligonucleotide sequences derived from well-known
FISH (fluorescence in situ hybridisation) probes and conserved
16S rRNA regions. Classification is hierarchical:

  1. Domain: Bacteria vs Archaea (EUB338 / archaeal markers)
  2. Phylum: within Bacteria (LGC354, HGC236, etc.)
  3. Class:  within Proteobacteria

Each probe entry stores the TARGET SITE in the 16S gene (sense strand),
i.e. the reverse complement of the published probe oligonucleotide.

IMPORTANT: Proteobacterial classes share nearly identical 16S V6 regions.
The key diagnostic nucleotide at position ~971 (ALF968 binding site) is:
  - Alpha: C  (AACGCG-C-AGAACCTTACC)
  - Beta/Gamma: G  (AACGCG-A/G-AGAACCTTACC)

Discrimination therefore requires EXACT matching (0 mismatches) for
class-specific probes, plus additional variable-region diagnostics.
"""

from __future__ import annotations

# ────────────────────────────────────────────────
# Domain-level probes
# ────────────────────────────────────────────────
# EUB338 target site (Amann et al. 1990) — present in >90% of Bacteria
# Probe: 5'-GCTGCCTCCCGTAGGAGT-3'  →  target in 16S gene:
BACTERIAL_MARKERS = [
    "ACTCCTACGGGAGGCAGC",   # EUB338 target (~pos 338)
    "ACTCCTACGGGAGGCAGC"[:12],  # core 12-mer for fuzzy matching
]

# Archaeal markers — sequences characteristic of archaeal 16S
# These are found in archaeal 16S but not bacterial 16S
ARCHAEAL_MARKERS = [
    "TCCGGTTGATCC",         # Archaeal 21F primer target region
    "ATTGGCGGGGG",          # ARC915 target core region
    "CCATGCGAGC",           # Archaeal-specific V3 region motif
]

# ────────────────────────────────────────────────
# Phylum / Class-level probes (for Bacteria only)
# ────────────────────────────────────────────────
#
# Proteobacterial class probes use max_mismatches=0 because
# the diagnostic nucleotides differ by only 1-2 bases between
# Alpha/Beta/Gamma. Even 1 mismatch causes cross-matching.
#
TAXONOMY_PROBES: dict[str, dict] = {
    "Firmicutes": {
        "description": "Firmicutes / Bacillota (Bacillus, Staphylococcus, Clostridium)",
        "domain": "Bacteria",
        "level": "phylum",
        # LGC354 a/b/c target sites (Meier et al. 1999) — specific to Firmicutes
        "probes": [
            "GCAGTAGGGAATCTTCCA",   # LGC354a target
            "GCAGTAGGGAATCTTCCT",   # LGC354b target
            "GCAGTAGGGAATCTCCCA",   # LGC354c target
        ],
        "max_mismatches": 2,
    },
    "Actinobacteria": {
        "description": "Actinobacteria (Streptomyces, Mycobacterium, Corynebacterium)",
        "domain": "Bacteria",
        "level": "phylum",
        # HGC236 target site — specific to high-GC Gram-positives
        "probes": [
            "TAACCGTATCAGCTTGTT",   # HGC236 target
        ],
        "max_mismatches": 2,
    },
    "Alphaproteobacteria": {
        "description": "Alpha-Proteobacteria (Rhizobium, Caulobacter, Rickettsia)",
        "domain": "Bacteria",
        "level": "class",
        # ALF968 target site — the C at position 7 is Alpha-diagnostic
        # (Beta/Gamma have G at this position)
        # Must be EXACT match — even 1 mismatch causes cross-matching with Gamma
        "probes": [
            "AACGCGCAGAACCTTACC",   # ALF968 target — C at pos 7 is key
        ],
        "max_mismatches": 0,  # MUST be exact — 1mm matches E. coli!
    },
    "Betaproteobacteria": {
        "description": "Beta-Proteobacteria (Burkholderia, Neisseria, Ralstonia)",
        "domain": "Bacteria",
        "level": "class",
        # Beta-specific V3 region — distinct from Gamma
        "probes": [
            "CCGCATACGCCCTTTGTAC",  # Beta-specific V3 helix
        ],
        "max_mismatches": 1,
    },
    "Gammaproteobacteria": {
        "description": "Gamma-Proteobacteria (E. coli, Pseudomonas, Vibrio, Salmonella)",
        "domain": "Bacteria",
        "level": "class",
        # Gamma-specific probes validated against E. coli K-12 MG1655 16S:
        #   pos 964: AACGCGAAGAACCTTACC (V6, shared with Beta — use exact)
        #   pos 767: GTGGGGAGCAAACAGG (ENT183 region, verified E. coli)
        # The V6 probe at pos 964 has G where Alpha has C — diagnostic
        "probes": [
            "AACGCGAAGAACCTTACC",   # V6 pos 964 — G at pos 7 (Alpha has C)
            "GTGGGGAGCAAACAGG",     # ENT183 region (corrected to actual E. coli)
            "CCTTTGTTGCCAGCG",      # Gamma V3 helix (validated match)
        ],
        "max_mismatches": 0,  # exact — V6 differs from Alpha by 1 base
    },
    "Deltaproteobacteria": {
        "description": "Delta-Proteobacteria (Myxococcus, Desulfovibrio, Bdellovibrio)",
        "domain": "Bacteria",
        "level": "class",
        "probes": [
            "CGGCGTCGCTGCGTCAGG",  # SRB385 target (sulfate reducers)
        ],
        "max_mismatches": 0,
    },
    "Bacteroidetes": {
        "description": "Bacteroidetes (Bacteroides, Flavobacterium, Cytophaga)",
        "domain": "Bacteria",
        "level": "phylum",
        # CF319 target site (Manz et al. 1996)
        "probes": [
            "TGGTCCGTGTCTCAGTAC",   # CF319a target
        ],
        "max_mismatches": 2,
    },
    "Cyanobacteria": {
        "description": "Cyanobacteria (Synechococcus, Nostoc, Anabaena)",
        "domain": "Bacteria",
        "level": "phylum",
        # CYA361 target site
        "probes": [
            "CCCATTGCGGAAAATTCC",   # CYA361 target
        ],
        "max_mismatches": 2,
    },
    "Spirochaetes": {
        "description": "Spirochaetes (Borrelia, Treponema, Leptospira)",
        "domain": "Bacteria",
        "level": "phylum",
        "probes": [
            "GCTTCGCAGGTGACTTTC",   # Spirochaete-diagnostic V3
        ],
        "max_mismatches": 2,
    },
}

# Archaeal subgroups (for classification within Archaea)
ARCHAEAL_PROBES: dict[str, dict] = {
    "Euryarchaeota": {
        "description": "Euryarchaeota (Methanogens, Halobacterium, Thermoplasma)",
        "level": "phylum",
        # EURY498 target site
        "probes": [
            "CTTGCCCRGCCCTT",       # EURY498 target (R = A or G)
            "CTTGCCCCGCCCTT",       # EURY498 variant A
            "CTTGCCCAGCCCTT",       # EURY498 variant G
        ],
        "max_mismatches": 2,
    },
    "Crenarchaeota": {
        "description": "Crenarchaeota (Sulfolobus, Thermoproteus)",
        "level": "phylum",
        # CREN499 target site
        "probes": [
            "CCAGRCTTGCCCCCC",      # CREN499 target
            "CCAGACTTGCCCCCC",      # variant A
            "CCAGGCTTGCCCCCC",      # variant G
        ],
        "max_mismatches": 2,
    },
}
