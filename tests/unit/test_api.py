"""Tests for the FastAPI REST backend."""

from __future__ import annotations

import pytest
from httpx import ASGITransport, AsyncClient

from darwin.api.app import app


@pytest.fixture
async def client():
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as ac:
        yield ac


class TestHealthEndpoint:
    async def test_health(self, client: AsyncClient):
        resp = await client.get("/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "ok"
        assert "version" in data
        assert "tools" in data


class TestAnnotateEndpoint:
    async def test_reject_non_fasta(self, client: AsyncClient):
        resp = await client.post(
            "/annotate",
            files={"file": ("test.txt", b"not a fasta", "text/plain")},
        )
        assert resp.status_code == 400

    async def test_submit_fasta(self, client: AsyncClient, tmp_fasta):
        with open(tmp_fasta, "rb") as f:
            resp = await client.post(
                "/annotate",
                files={"file": (tmp_fasta.name, f, "application/octet-stream")},
            )
        assert resp.status_code == 200
        data = resp.json()
        assert "job_id" in data
        assert data["status"] == "submitted"


class TestJobsEndpoint:
    async def test_job_not_found(self, client: AsyncClient):
        resp = await client.get("/jobs/nonexistent")
        assert resp.status_code == 404
