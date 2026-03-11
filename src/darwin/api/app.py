"""Darwin REST API — FastAPI backend for genome annotation.

Endpoints:
    POST /annotate          — Upload FASTA, get annotation results
    GET  /jobs/{job_id}     — Check job status / retrieve results
    GET  /health            — Health check
    GET  /docs              — Auto-generated Swagger UI
"""

from __future__ import annotations

import shutil
import tempfile
import uuid
from pathlib import Path
from typing import Any

from fastapi import BackgroundTasks, FastAPI, File, HTTPException, Query, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel

from darwin import __version__
from darwin.api.jobs import JobManager, JobStatus

app = FastAPI(
    title="Darwin Genome Annotator",
    description="Fast, accurate prokaryotic genome annotation API",
    version=__version__,
)

# CORS — allow frontend from any origin (tighten in production)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# In-memory job manager (swap for Redis/DB in production)
jobs = JobManager()


# ── Request / Response Models ────────────────────────────────


class AnnotateParams(BaseModel):
    locus_tag: str = "DARWIN"
    translation_table: int = 11
    evalue: float = 1e-6
    metagenome: bool = False
    kingdom: str = "bac"
    min_contig_len: int = 200
    cpus: int = 1


class JobResponse(BaseModel):
    job_id: str
    status: str
    message: str


class JobResultResponse(BaseModel):
    job_id: str
    status: str
    result: dict[str, Any] | None = None
    error: str | None = None


class HealthResponse(BaseModel):
    status: str
    version: str
    tools: dict[str, bool]


# ── Endpoints ────────────────────────────────────────────────


@app.get("/health", response_model=HealthResponse)
async def health_check() -> HealthResponse:
    """Check API health and tool availability."""
    from darwin.annotators import (
        AragornAnnotator,
        BarrnapAnnotator,
        ProdigalAnnotator,
        PyhmmerAnnotator,
    )
    from darwin.models import AnnotationConfig

    dummy = AnnotationConfig(
        input_path=Path("/dev/null"),
        output_dir=Path(tempfile.mkdtemp()),
    )

    tools = {
        "prodigal": ProdigalAnnotator(dummy).check_dependencies(),
        "barrnap": BarrnapAnnotator(dummy).check_dependencies(),
        "aragorn": AragornAnnotator(dummy).check_dependencies(),
        "pyhmmer": PyhmmerAnnotator(dummy).check_dependencies(),
    }

    return HealthResponse(
        status="ok",
        version=__version__,
        tools=tools,
    )


@app.post("/annotate", response_model=JobResponse)
async def annotate_genome(
    background_tasks: BackgroundTasks,
    file: UploadFile = File(..., description="FASTA file to annotate"),
    locus_tag: str = Query("DARWIN", description="Locus tag prefix"),
    translation_table: int = Query(11, description="NCBI translation table"),
    evalue: float = Query(1e-6, description="HMM e-value threshold"),
    metagenome: bool = Query(False, description="Metagenomic mode"),
    kingdom: str = Query("bac", description="Kingdom: bac or arc"),
    min_contig_len: int = Query(200, description="Min contig length"),
    cpus: int = Query(1, description="CPU threads"),
) -> JobResponse:
    """Submit a genome for annotation.

    Returns a job_id to poll for results.
    """
    if not file.filename or not file.filename.endswith((".fasta", ".fa", ".fna", ".fasta.gz")):
        raise HTTPException(
            status_code=400,
            detail="File must be FASTA format (.fasta, .fa, .fna)",
        )

    # Create job
    job_id = str(uuid.uuid4())[:8]
    work_dir = Path(tempfile.mkdtemp(prefix=f"darwin_{job_id}_"))
    input_path = work_dir / file.filename

    # Save uploaded file
    content = await file.read()
    input_path.write_bytes(content)

    # Register job
    jobs.create(job_id)

    # Run annotation in background
    from darwin.models import AnnotationConfig

    config = AnnotationConfig(
        input_path=input_path,
        output_dir=work_dir / "output",
        locus_tag_prefix=locus_tag,
        translation_table=translation_table,
        evalue=evalue,
        cpus=cpus,
        metagenome=metagenome,
        kingdom=kingdom,
        min_contig_len=min_contig_len,
        formats=["json"],  # API only needs JSON
    )

    background_tasks.add_task(_run_annotation, job_id, config)

    return JobResponse(
        job_id=job_id,
        status="submitted",
        message=f"Annotation job {job_id} submitted. Poll GET /jobs/{job_id} for results.",
    )


@app.get("/jobs/{job_id}", response_model=JobResultResponse)
async def get_job(job_id: str) -> JobResultResponse:
    """Get the status and results of an annotation job."""
    job = jobs.get(job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    return JobResultResponse(
        job_id=job_id,
        status=job["status"],
        result=job.get("result"),
        error=job.get("error"),
    )


# ── Background Task ──────────────────────────────────────────


def _run_annotation(job_id: str, config: AnnotationConfig) -> None:
    """Run the annotation pipeline as a background task."""
    from darwin.output.json_out import genome_to_dict
    from darwin.pipeline.runner import DarwinPipeline

    jobs.update(job_id, status="running")

    try:
        pipeline = DarwinPipeline(config)
        genome = pipeline.run()
        result = genome_to_dict(genome)
        jobs.update(job_id, status="completed", result=result)
    except Exception as exc:
        jobs.update(job_id, status="failed", error=str(exc))
    finally:
        # Clean up temp files (keep results in memory)
        try:
            shutil.rmtree(config.output_dir.parent, ignore_errors=True)
        except Exception:
            pass
