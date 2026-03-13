"""
Downloader — Automated database retrieval for the soil.

Like rain falling into the jar — it brings essential nutrients
(HMM databases) from the outside world so organisms can thrive.
"""

from __future__ import annotations

import gzip
import logging
import shutil
import subprocess
import urllib.request
from pathlib import Path

logger = logging.getLogger("darwin.soil.downloader")

# Default cache location
DEFAULT_CACHE_DIR = Path.home() / ".darwin" / "databases"

# Database URLs — with fallback mirrors for reliability
DATABASES = {
    "pfam": {
        "name": "Pfam-A",
        "urls": [
            "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
            "http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
        ],
        "filename": "Pfam-A.hmm",
        "description": "Pfam-A protein families (EBI/InterPro)",
    },
    "tigrfams": {
        "name": "TIGRFAMs",
        "urls": [
            "https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release/TIGRFAMs_15.0_HMM.LIB.gz",
            "https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.LIB.gz",
        ],
        "filename": "TIGRFAMs.hmm",
        "description": "TIGRFAMs protein families (NCBI/JCVI)",
    },
}


def get_cache_dir() -> Path:
    """Get the database cache directory, creating it if needed."""
    DEFAULT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return DEFAULT_CACHE_DIR


def list_cached_databases() -> list[Path]:
    """Find all HMM databases in the cache directory."""
    cache = get_cache_dir()
    dbs = []
    for hmm_file in sorted(cache.glob("*.hmm")):
        dbs.append(hmm_file)
    return dbs


def download_database(
    db_key: str,
    cache_dir: Path | None = None,
    force: bool = False,
) -> Path | None:
    """
    Download and prepare an HMM database.

    Returns the path to the ready-to-use .hmm file, or None on failure.
    """
    if db_key not in DATABASES:
        logger.error(f"Unknown database: {db_key}. Available: {list(DATABASES.keys())}")
        return None

    db_info = DATABASES[db_key]
    target_dir = cache_dir or get_cache_dir()
    target_dir.mkdir(parents=True, exist_ok=True)

    hmm_path = target_dir / db_info["filename"]

    # Check if already exists and pressed
    if hmm_path.exists() and hmm_path.with_suffix(".hmm.h3i").exists() and not force:
        logger.info(f"✅ {db_info['name']} already cached at {hmm_path}")
        return hmm_path

    # Download — try each URL in order until one succeeds
    urls = db_info.get("urls", [db_info.get("url", "")])
    gz_path = hmm_path.with_suffix(".hmm.gz")

    downloaded = False
    for url in urls:
        logger.info(f"⬇️  Downloading {db_info['name']} from {url}...")
        try:
            _download_with_progress(url, gz_path)
            downloaded = True
            break
        except Exception as e:
            logger.warning(f"Mirror failed ({e}), trying next...")

    if not downloaded:
        logger.error(
            f"All download mirrors failed for {db_info['name']}. "
            f"You can manually download and place the .hmm file at: {hmm_path}"
        )
        return None

    # Decompress
    logger.info(f"📦 Decompressing {gz_path.name}...")
    try:
        with gzip.open(gz_path, "rb") as f_in, open(hmm_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        gz_path.unlink()  # Remove compressed file
    except Exception as e:
        logger.error(f"Decompression failed: {e}")
        return None

    # Press with hmmpress
    logger.info(f"🔨 Pressing {hmm_path.name} with hmmpress...")
    hmmpress = shutil.which("hmmpress")
    if not hmmpress:
        logger.warning(
            "hmmpress not found — database downloaded but not indexed. "
            "Install HMMER and run: hmmpress " + str(hmm_path)
        )
        return hmm_path

    try:
        result = subprocess.run(
            [hmmpress, "-f", str(hmm_path)],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode != 0:
            logger.error(f"hmmpress failed: {result.stderr}")
            return hmm_path  # Still usable by pyhmmer without pressing
    except subprocess.TimeoutExpired:
        logger.error("hmmpress timed out")
        return hmm_path

    logger.info(f"✅ {db_info['name']} ready at {hmm_path}")
    return hmm_path


def download_all(
    cache_dir: Path | None = None,
    force: bool = False,
) -> list[Path]:
    """Download all available databases."""
    paths = []
    for key in DATABASES:
        path = download_database(key, cache_dir=cache_dir, force=force)
        if path:
            paths.append(path)
    return paths


def _download_with_progress(url: str, dest: Path) -> None:
    """Download a file with basic progress reporting."""
    req = urllib.request.Request(url, headers={"User-Agent": "Darwin-Annotator/0.2"})

    with urllib.request.urlopen(req, timeout=600) as response:
        total = int(response.headers.get("Content-Length", 0))
        downloaded = 0
        chunk_size = 1024 * 1024  # 1MB chunks

        with open(dest, "wb") as f:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)

                if total > 0:
                    pct = downloaded / total * 100
                    mb_done = downloaded / 1024 / 1024
                    mb_total = total / 1024 / 1024
                    logger.info(f"  {mb_done:.0f}/{mb_total:.0f} MB ({pct:.0f}%)")

    logger.info(f"  Downloaded {downloaded / 1024 / 1024:.1f} MB")
