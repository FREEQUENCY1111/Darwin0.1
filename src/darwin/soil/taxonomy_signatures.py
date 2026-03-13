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
TAXONOMY_PROBES: dict[str, dict] = {
    "Firmicutes": {
        "description": "Firmicutes / Bacillota (Bacillus, Staphylococcus, Clostridium)",
        "domain": "Bacteria",
        "level": "phylum",
        # LGC354 a/b/c target sites (Meier et al. 1999) — specific to Firmicutes
        # Probe LGC354a: 5'-TGGAAGATTCCCTACTGC-3'  → target:
        "probes": [
            "GCAGTAGGGAATCTTCCA",   # LGC354a target
            "GCAGTAGGGAATCTTCCT",   # LGC354b target
            "GCAGTAGGGAATCTCCCA",   # LGC354c target
        ],
        "min_hits": 1,
    },
    "Actinobacteria": {
        "description": "Actinobacteria (Streptomyces, Mycobacterium, Corynebacterium)",
        "domain": "Bacteria",
        "level": "phylum",
        # HGC236 target site — specific to high-GC Gram-positives
        # Probe HGC236: 5'-AACAAGCTGATACGGTTA-3'  → target:
        "probes": [
            "TAACCGTATCAGCTTGTT",   # HGC236 target
        ],
        "min_hits": 1,
    },
    "Alphaproteobacteria": {
        "description": "Alpha-Proteobacteria (Rhizobium, Caulobacter, Rickettsia)",
        "domain": "Bacteria",
        "level": "class",
        # ALF968 target site (Neef et al. 1999)
        # Probe ALF968: 5'-GGTAAGGTTCTGCGCGTT-3'  → target:
        "probes": [
            "AACGCGCAGAACCTTACC",   # ALF968 target
        ],
        "min_hits": 1,
    },
    "Betaproteobacteria": {
        "description": "Beta-Proteobacteria (Burkholderia, Neisseria, Ralstonia)",
        "domain": "Bacteria",
        "level": "class",
        # BET42a is a 23S probe, so use V6 diagnostic motifs
        "probes": [
            "AACGCGAAGAACCTTACC",   # Beta-specific V6 variant
        ],
        "min_hits": 1,
    },
    "Gammaproteobacteria": {
        "description": "Gamma-Proteobacteria (E. coli, Pseudomonas, Vibrio)",
        "domain": "Bacteria",
        "level": "class",
        # Gammaproteobacteria have distinctive V3 and V6 regions
        # GAM42a is 23S, so use 16S diagnostic motifs
        "probes": [
            "AACGCGAAAAACCTTACC",   # Gamma V6 variant
            "AACGCGAAGAACCTTACCT",  # Gamma extended
        ],
        "min_hits": 1,
    },
    "Bacteroidetes": {
        "description": "Bacteroidetes (Bacteroides, Flavobacterium, Cytophaga)",
        "domain": "Bacteria",
        "level": "phylum",
        # CF319 target site (Manz et al. 1996)
        "probes": [
            "TGGTCCGTGTCTCAGTAC",   # CF319a target
        ],
        "min_hits": 1,
    },
    "Cyanobacteria": {
        "description": "Cyanobacteria (Synechococcus, Nostoc, Anabaena)",
        "domain": "Bacteria",
        "level": "phylum",
        # CYA361 target site
        "probes": [
            "CCCATTGCGGAAAATTCC",   # CYA361 target
        ],
        "min_hits": 1,
    },
    "Spirochaetes": {
        "description": "Spirochaetes (Borrelia, Treponema, Leptospira)",
        "domain": "Bacteria",
        "level": "phylum",
        "probes": [
            "GCTTCGCAGGTGACTTTC",   # Spirochaete-diagnostic V3
        ],
        "min_hits": 1,
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
        "min_hits": 1,
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
        "min_hits": 1,
    },
}
