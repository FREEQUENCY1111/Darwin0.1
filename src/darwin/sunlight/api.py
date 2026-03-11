"""
API — Artificial sunlight. Same energy, different source.

FastAPI endpoints that let you interact with the ecosphere
over HTTP. Like growing plants under grow lights instead
of natural sun — same result, different delivery.
"""

from __future__ import annotations

import asyncio
import time
import uuid
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException, UploadFile, File, Form
from pydantic import BaseModel

from darwin.rocks.models import AnnotationConfig
from darwin.jar.ecosphere import Ecosphere


class JobStatus(BaseModel):
    job_id: str
    status: str  # pending, running, completed, failed
    created_at: float
    completed_at: Optional[float] = None
    result: Optional[dict] = None
    error: Optional[str] = None


# In-memory job store (swap for Redis in production)
_jobs: dict[str, JobStatus] = {}


def create_app() -> FastAPI:
    """Create the Darwin API — artificial sunlight."""

    app = FastAPI(
        title="Darwin Ecosphere API",
        description=(
            "Prokaryotic genome annotation powered by a self-sustaining "
            "ecosphere. Just add sunlight (your FASTA file)."
        ),
        version="0.1.0",
    )

    @app.get("/")
    async def root():
        return {
            "name": "Darwin Ecosphere",
            "status": "alive",
            "message": "Add sunlight to begin. POST /annotate with a FASTA file.",
        }

    @app.get("/health")
    async def health():
        """Check ecosystem health — is there enough soil?"""
        from darwin.soil.nutrients import NutrientStore
        soil = NutrientStore()
        report = soil.survey()
        return {
            "status": "healthy" if soil.is_fertile else "barren",
            "soil": report,
        }

    @app.post("/annotate")
    async def annotate(
        file: UploadFile = File(...),
        output_dir: str = Form("darwin_output"),
        locus_prefix: str = Form("DARWIN"),
        translation_table: int = Form(11),
        evalue: float = Form(1e-10),
        metagenome: bool = Form(False),
    ):
        """
        Add sunlight — submit a genome for annotation.

        Returns a job ID. Poll /jobs/{job_id} for results.
        The ecosphere works asynchronously — equilibrium
        takes time.
        """
        job_id = str(uuid.uuid4())[:8]

        # Save uploaded file
        upload_dir = Path("/tmp/darwin_uploads")
        upload_dir.mkdir(exist_ok=True)
        fasta_path = upload_dir / f"{job_id}_{file.filename}"

        content = await file.read()
        with open(fasta_path, "wb") as fh:
            fh.write(content)

        # Create job
        _jobs[job_id] = JobStatus(
            job_id=job_id,
            status="pending",
            created_at=time.time(),
        )

        # Run ecosphere in background
        config = AnnotationConfig(
            input_file=fasta_path,
            output_dir=Path(output_dir) / job_id,
            locus_tag_prefix=locus_prefix,
            translation_table=translation_table,
            evalue_threshold=evalue,
            metagenome_mode=metagenome,
        )

        asyncio.create_task(_run_ecosphere(job_id, config, fasta_path))

        return {
            "job_id": job_id,
            "status": "pending",
            "message": "☀️ Sunlight added. Ecosystem is growing...",
            "poll_url": f"/jobs/{job_id}",
        }

    @app.get("/jobs/{job_id}")
    async def get_job(job_id: str):
        """Check on the ecosystem — has equilibrium been reached?"""
        job = _jobs.get(job_id)
        if not job:
            raise HTTPException(404, f"No ecosystem found with id {job_id}")
        return job

    @app.get("/jobs")
    async def list_jobs():
        """List all ecosystems (past and present)."""
        return {
            "jobs": [
                {"job_id": j.job_id, "status": j.status, "created_at": j.created_at}
                for j in _jobs.values()
            ]
        }

    return app


async def _run_ecosphere(
    job_id: str, config: AnnotationConfig, fasta_path: Path
) -> None:
    """Run the ecosphere in the background."""
    _jobs[job_id].status = "running"

    try:
        jar = Ecosphere(config)
        result = await jar.add_sunlight(fasta_path=fasta_path)

        _jobs[job_id].status = "completed"
        _jobs[job_id].completed_at = time.time()
        _jobs[job_id].result = result

    except Exception as e:
        _jobs[job_id].status = "failed"
        _jobs[job_id].completed_at = time.time()
        _jobs[job_id].error = str(e)

    finally:
        # Clean up uploaded file
        try:
            fasta_path.unlink(missing_ok=True)
        except Exception:
            pass
