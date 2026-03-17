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

Multiple probes per group improve discrimination — especially for
Proteobacterial classes, which share very similar V6 regions and
cannot be reliably separated by a single probe.
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
# For Proteobacterial classes, 16S probes alone are insufficient
# because the V6 region (around pos 968) is nearly identical
# across Alpha/Beta/Gamma. We therefore use probes from multiple
# variable regions (V1, V3, V6) and require stricter matching
# (max 1 mismatch) for proteobacterial class-level probes.
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
        # ALF968 target site (Neef et al. 1999) plus Alpha-specific V1 region
        "probes": [
            "AACGCGCAGAACCTTACC",   # ALF968 target — note C at pos 7
            "CGGTAATACGGAGGGTGC",   # Alpha V1 signature
            "GTAGTCCACGCCGTAAAC",   # Alpha V4 diagnostic
        ],
        "max_mismatches": 1,  # strict — proteobacterial classes are similar
    },
    "Betaproteobacteria": {
        "description": "Beta-Proteobacteria (Burkholderia, Neisseria, Ralstonia)",
        "domain": "Bacteria",
        "level": "class",
        # Beta-specific V6 and V3 regions
        "probes": [
            "AACGCGAAGAACCTTACC",   # Beta V6 — note shared with Gamma
            "CCGCATACGCCCTTTGTAC",  # Beta-specific V3 helix
            "GGAATCTTGCGCAATGGG",  # Beta V1 diagnostic
        ],
        "max_mismatches": 1,
    },
    "Gammaproteobacteria": {
        "description": "Gamma-Proteobacteria (E. coli, Pseudomonas, Vibrio, Salmonella)",
        "domain": "Bacteria",
        "level": "class",
        # Gamma-specific probes from multiple 16S regions
        # The V6 region is shared with Beta, so discrimination
        # relies on V3 (Enterobacteriaceae/Pseudomonadales) and V1 signatures
        "probes": [
            "GTGGGGAGCAAAGAGC",     # ENT183 target — Enterobacteriaceae V3
            "CCTTTGTTGCCAGCG",      # Gamma V3 helix signature
            "ATGACGGTACCTGAGAAG",   # Gamma V4 diagnostic
            "GGGAGTACGGTCGCAAG",    # Gamma V1 signature (broad)
        ],
        "max_mismatches": 1,
    },
    "Deltaproteobacteria": {
        "description": "Delta-Proteobacteria (Myxococcus, Desulfovibrio, Bdellovibrio)",
        "domain": "Bacteria",
        "level": "class",
        "probes": [
            "CGGCGTCGCTGCGTCAGG",  # SRB385 target (sulfate reducers)
        ],
        "max_mismatches": 1,
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
