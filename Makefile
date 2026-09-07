.PHONY: help clean clean-build clean-pyc clean-test lint format typecheck test test-cov \
        coverage build dist preflight docs servedocs install \
        bump-patch bump-minor bump-major
.DEFAULT_GOAL := help

# Every target runs through one interpreter. Override for a venv or a specific
# version:  make test PYTHON=.venv/bin/python
#
# Tools are invoked as `$(PYTHON) -m <tool>`, never as bare executables: a
# `pip install --user` puts them in a scripts dir that is usually not on PATH,
# which is how `make bump-*` fails with "command not found" on a machine where
# the tool is definitely installed. The module form is equivalent and
# PATH-independent.
PYTHON ?= python3

help:
	@echo "PYTHON=$(PYTHON)"
	@echo ""
	@echo "clean        remove build, test and Python artifacts"
	@echo "lint         run ruff + black --check (the CI lint gate)"
	@echo "format       apply black + ruff --fix"
	@echo "typecheck    run mypy (fatal, as in CI)"
	@echo "test         run the test suite"
	@echo "test-cov     run tests with coverage"
	@echo "build/dist   build sdist + wheel"
	@echo "preflight    build + twine check + wheel-name check (what CI asserts)"
	@echo "docs         build the Sphinx HTML docs"
	@echo "bump-patch   bump patch version, commit and tag (vX.Y.Z)"
	@echo "bump-minor   bump minor version, commit and tag"
	@echo "bump-major   bump major version, commit and tag"

clean: clean-build clean-pyc clean-test

clean-build:
	rm -fr build/ dist/ .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc:
	find . -name '*.pyc' -delete
	find . -name '*.pyo' -delete
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -fr .tox/ .pytest_cache .mypy_cache htmlcov/
	rm -f .coverage coverage.xml

# ruff replaces flake8 *and* isort (REPO_CONTRACT §4b). Do not add either back.
lint:
	$(PYTHON) -m ruff check RiboMetric tests
	$(PYTHON) -m black --check RiboMetric tests

format:
	$(PYTHON) -m black RiboMetric tests
	$(PYTHON) -m ruff check --fix RiboMetric tests

typecheck:
	$(PYTHON) -m mypy --config-file mypy.ini RiboMetric

test:
	RIBOMETRIC_SKIP_IMAGES=1 $(PYTHON) -m pytest -q -n auto

test-cov coverage:
	RIBOMETRIC_SKIP_IMAGES=1 $(PYTHON) -m pytest -n auto \
		--cov=RiboMetric --cov-report=term-missing --cov-report=xml

build dist: clean
	$(PYTHON) -m build
	ls -l dist

# Mirrors ci.yml's `build` job exactly, and is the only rehearsal for a
# release. There is deliberately no TestPyPI target (REPO_CONTRACT §5.4): it
# uses a stored token and exercises neither OIDC nor the trusted-publisher
# binding, so it cannot rehearse the real path. `twine check` alone does not
# catch a non-normalised wheel filename -- the check below is what does.
preflight: dist
	$(PYTHON) -m twine check dist/*
	@for f in dist/*.whl dist/*.tar.gz; do \
		case "$$(basename $$f)" in \
			ribometric-*) ;; \
			*) echo "ERROR: $$(basename $$f) is not normalised; PyPI needs the 'ribometric-' prefix"; exit 1 ;; \
		esac; \
	done
	@echo "dist filenames normalised:"; ls -1 dist

docs:
	rm -f docs/RiboMetric.rst docs/modules.rst
	$(PYTHON) -m sphinx.ext.apidoc -o docs/ RiboMetric
	$(MAKE) -C docs clean
	$(MAKE) -C docs html

servedocs: docs
	$(PYTHON) -m http.server --directory docs/_build/html

install: clean
	$(PYTHON) -m pip install -e '.[dev]'

# bump2version is the only way the version changes (REPO_CONTRACT §5.2). It
# rewrites every file in .bumpversion.cfg, commits and tags in one step, so the
# bump *is* the release-cut action. There is no `release` target that uploads:
# pushing the tag is what publishes, through release.yml.
bump-patch:
	$(PYTHON) -m bumpversion patch

bump-minor:
	$(PYTHON) -m bumpversion minor

bump-major:
	$(PYTHON) -m bumpversion major
