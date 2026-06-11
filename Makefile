.PHONY: help install install-dev test test-unit test-integration test-benchmark coverage lint format clean docs

help:
	@echo "Available commands:"
	@echo "  make install          Install package"
	@echo "  make install-dev      Install package with dev dependencies"
	@echo "  make test             Run all tests"
	@echo "  make test-unit        Run unit tests only"
	@echo "  make test-integration Run integration tests only"
	@echo "  make test-benchmark   Run benchmark tests only"
	@echo "  make coverage         Generate coverage report"
	@echo "  make lint             Run linters (flake8, mypy)"
	@echo "  make format           Format code (black, isort)"
	@echo "  make clean            Clean build artifacts"
	@echo "  make docs             Build documentation"

install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"

test:
	pytest tests/

test-unit:
	pytest tests/unit/ -m unit

test-integration:
	pytest tests/integration/ -m integration

test-benchmark:
	pytest tests/benchmark/ -m benchmark

coverage:
	pytest --cov=spagapa --cov-report=html --cov-report=term
	@echo "Coverage report generated in htmlcov/index.html"

lint:
	flake8 spagapa tests
	mypy spagapa

format:
	black spagapa tests
	isort spagapa tests

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache
	rm -rf .coverage
	rm -rf htmlcov/
	rm -rf .mypy_cache
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

docs:
	cd docs && make html
	@echo "Documentation built in docs/_build/html/index.html"
