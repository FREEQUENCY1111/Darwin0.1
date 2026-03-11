.PHONY: install dev test lint typecheck check serve docker clean

# ── Development ──────────────────────────────────────────────

install:
	pip install -e ".[dev,test]"

dev: install
	pre-commit install
	darwin check

# ── Quality ──────────────────────────────────────────────────

lint:
	ruff check src/ tests/
	ruff format --check src/ tests/

format:
	ruff check --fix src/ tests/
	ruff format src/ tests/

typecheck:
	mypy src/ --ignore-missing-imports

test:
	pytest tests/unit/ -v --tb=short --cov=src --cov-report=term-missing

test-all:
	pytest tests/ -v --tb=short --cov=src --cov-report=xml

# ── Running ──────────────────────────────────────────────────

check:
	darwin check

serve:
	darwin serve --reload --port 8000

# ── Docker ───────────────────────────────────────────────────

docker:
	docker compose up --build

docker-bg:
	docker compose up --build -d

docker-down:
	docker compose down

# ── Cleanup ──────────────────────────────────────────────────

clean:
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .pytest_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .mypy_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .ruff_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	rm -rf dist/ build/ coverage.xml htmlcov/
