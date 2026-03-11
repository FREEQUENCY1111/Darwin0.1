.PHONY: install dev lint test serve docker

install:
	pip install -e .

dev:
	pip install -e ".[dev,test]"

lint:
	ruff check src/ tests/
	ruff format --check src/ tests/

format:
	ruff format src/ tests/

test:
	pytest tests/ -v --tb=short

test-cov:
	pytest tests/ -v --cov=darwin --cov-report=html

serve:
	darwin serve --port 8000

check:
	darwin check

docker:
	docker compose up --build
