# ============================================================
# Dockerfile — Prokaryotic Genome Annotator
# ============================================================
# Multi-stage build for a slim, reproducible image with all
# bioinformatics dependencies pre-installed.

# ── Stage 1: Build ───────────────────────────────────────────
FROM python:3.11-slim AS builder

WORKDIR /build
COPY pyproject.toml README.md ./
COPY src/ src/

RUN pip install --no-cache-dir build \
    && python -m build --wheel

# ── Stage 2: Runtime ─────────────────────────────────────────
FROM python:3.11-slim

LABEL org.opencontainers.image.source="https://github.com/YOUR_USERNAME/genome-annotator"
LABEL org.opencontainers.image.description="Prokaryotic genome annotation pipeline"

# Install bioinformatics system dependencies
RUN apt-get update -q \
    && apt-get install -y --no-install-recommends \
       prodigal \
       hmmer \
       aragorn \
       barrnap \
       infernal \
    && rm -rf /var/lib/apt/lists/*

# Install the Python package from the built wheel
COPY --from=builder /build/dist/*.whl /tmp/
RUN pip install --no-cache-dir /tmp/*.whl \
    && rm /tmp/*.whl

# Create a non-root user for safety
RUN useradd --create-home annotator
USER annotator
WORKDIR /home/annotator

ENTRYPOINT ["genome-annotator"]
CMD ["--help"]
