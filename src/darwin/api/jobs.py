"""Simple in-memory job manager for the Darwin API.

For production, swap this with Redis or a database-backed store.
"""

from __future__ import annotations

import threading
from typing import Any


class JobStatus:
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class JobManager:
    """Thread-safe in-memory job store."""

    def __init__(self) -> None:
        self._jobs: dict[str, dict[str, Any]] = {}
        self._lock = threading.Lock()

    def create(self, job_id: str) -> None:
        with self._lock:
            self._jobs[job_id] = {
                "status": JobStatus.PENDING,
                "result": None,
                "error": None,
            }

    def update(
        self,
        job_id: str,
        *,
        status: str | None = None,
        result: dict | None = None,
        error: str | None = None,
    ) -> None:
        with self._lock:
            if job_id not in self._jobs:
                return
            if status:
                self._jobs[job_id]["status"] = status
            if result is not None:
                self._jobs[job_id]["result"] = result
            if error is not None:
                self._jobs[job_id]["error"] = error

    def get(self, job_id: str) -> dict[str, Any] | None:
        with self._lock:
            return self._jobs.get(job_id)

    def list_jobs(self) -> dict[str, str]:
        """Return {job_id: status} for all jobs."""
        with self._lock:
            return {jid: j["status"] for jid, j in self._jobs.items()}
