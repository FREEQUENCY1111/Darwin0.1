FROM python:3.11-slim AS base

# Install bioinformatics tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    prodigal \
    aragorn \
    barrnap \
    hmmer \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user
RUN useradd -m -s /bin/bash darwin
WORKDIR /app

# Install Python package
COPY pyproject.toml README.md ./
COPY src/ src/
RUN pip install --no-cache-dir -e "."

USER darwin

ENTRYPOINT ["darwin"]
CMD ["--help"]
