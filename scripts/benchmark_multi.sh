#!/bin/bash
# ============================================================
# Darwin v0.2 Multi-Organism Benchmark
# Tests: B. subtilis (Gram+ bacterium) + M. jannaschii (archaeon)
# Run from the DARWIN project root directory
# ============================================================
set -e

echo "============================================================"
echo "  Darwin v0.2 Multi-Organism Benchmark"
echo "============================================================"

# ---- Bacillus subtilis str. 168 (GCF_000009045.1) ----
# Gram-positive model bacterium, ~4.2 Mb, 4,243 CDS
echo ""
echo ">>> [1/6] Downloading B. subtilis str. 168..."
mkdir -p bsub_benchmark
cd bsub_benchmark
if [ ! -f GCF_000009045.1_ASM904v1_genomic.fna ]; then
    curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz" | gunzip > GCF_000009045.1_ASM904v1_genomic.fna
    curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.gff.gz" | gunzip > genomic.gff
    echo "  Downloaded $(wc -c < GCF_000009045.1_ASM904v1_genomic.fna) bytes"
else
    echo "  Already downloaded, skipping."
fi
cd ..

echo ""
echo ">>> [2/6] Annotating B. subtilis with Darwin..."
darwin annotate \
    bsub_benchmark/GCF_000009045.1_ASM904v1_genomic.fna \
    -o bsub_output \
    -p BSUB

echo ""
echo ">>> [3/6] Benchmarking B. subtilis..."
python scripts/benchmark.py \
    --ncbi bsub_benchmark/genomic.gff \
    --darwin bsub_output/GCF_000009045.1_ASM904v1_genomic.gff3 \
    --output benchmark_bsub.json

# ---- Methanocaldococcus jannaschii DSM 2661 (GCF_000091665.1) ----
# First sequenced archaeon, ~1.7 Mb, 1,786 CDS
echo ""
echo ">>> [4/6] Downloading M. jannaschii DSM 2661..."
mkdir -p mjan_benchmark
cd mjan_benchmark
if [ ! -f GCF_000091665.1_ASM9166v1_genomic.fna ]; then
    curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/665/GCF_000091665.1_ASM9166v1/GCF_000091665.1_ASM9166v1_genomic.fna.gz" | gunzip > GCF_000091665.1_ASM9166v1_genomic.fna
    curl -L "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/665/GCF_000091665.1_ASM9166v1/GCF_000091665.1_ASM9166v1_genomic.gff.gz" | gunzip > genomic.gff
    echo "  Downloaded $(wc -c < GCF_000091665.1_ASM9166v1_genomic.fna) bytes"
else
    echo "  Already downloaded, skipping."
fi
cd ..

echo ""
echo ">>> [5/6] Annotating M. jannaschii with Darwin..."
darwin annotate \
    mjan_benchmark/GCF_000091665.1_ASM9166v1_genomic.fna \
    -o mjan_output \
    -p MJAN

echo ""
echo ">>> [6/6] Benchmarking M. jannaschii..."
python scripts/benchmark.py \
    --ncbi mjan_benchmark/genomic.gff \
    --darwin mjan_output/GCF_000091665.1_ASM9166v1_genomic.gff3 \
    --output benchmark_mjan.json

echo ""
echo "============================================================"
echo "  All benchmarks complete!"
echo "  Results saved to: benchmark_bsub.json, benchmark_mjan.json"
echo "============================================================"
