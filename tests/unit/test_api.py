"""Tests for the Sunlight API layer."""

import pytest
from httpx import AsyncClient, ASGITransport

from darwin.sunlight.api import create_app


@pytest.fixture
def app():
    return create_app()


class TestAPI:
    @pytest.mark.asyncio
    async def test_root(self, app):
        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as client:
            resp = await client.get("/")
            assert resp.status_code == 200
            data = resp.json()
            assert "Darwin" in data["name"]

    @pytest.mark.asyncio
    async def test_health(self, app):
        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as client:
            resp = await client.get("/health")
            assert resp.status_code == 200
            data = resp.json()
            assert "status" in data
            assert "soil" in data

    @pytest.mark.asyncio
    async def test_jobs_empty(self, app):
        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as client:
            resp = await client.get("/jobs")
            assert resp.status_code == 200

    @pytest.mark.asyncio
    async def test_job_not_found(self, app):
        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://test") as client:
            resp = await client.get("/jobs/nonexistent")
            assert resp.status_code == 404
