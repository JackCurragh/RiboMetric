#!/usr/bin/env python

"""The setup script."""

import re
from pathlib import Path

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

# Single source of truth for the version: RiboMetric/__init__.py
_init = Path(__file__).parent / "RiboMetric" / "__init__.py"
version = re.search(
    r'^__version__\s*=\s*["\']([^"\']+)["\']', _init.read_text(), re.M
).group(1)

# Core dependencies needed to run RiboMetric
requirements = [
    "biopython>=1.81",
    "gffpandas>=1.2.0",
    "Jinja2>=3.0",
    # Prefer newer Kaleido (>=1.0.0) to avoid deprecations and get wheels
    "kaleido>=1.0.0",  # Required for plotly image export
    # numpy<2 / pandas<2.3: upper caps reflect the tested-compatible range with
    # oxbow/pyarrow. Revisit once validated against numpy 2.x.
    "numpy>=1.24,<2",
    "oxbow>=0.3.0",
    "pandas>=2.0,<2.3",
    "plotly>=5.14",
    "pyarrow>=10.0.0",
    "pysam>=0.21.0",
    "PyYAML>=6.0",
    "rich>=13.3.3",  # Required by textual
    "scipy>=1.10",
    "textual>=0.47.0",
]

# Optional PDF export dependencies
pdf_requirements = [
    "xhtml2pdf>=0.2.11",
]

# Test dependencies
test_requirements = [
    "pytest>=7.0.0",
    "pytest-cov>=4.0.0",
    "pytest-mock>=3.10.0",
    "pytest-xdist>=3.0.0",
    "coverage>=7.0.0",
]

# Development dependencies (linting, formatting, docs)
# ruff replaces flake8 *and* isort (REPO_CONTRACT §4b); do not add either back.
dev_requirements = test_requirements + [
    "ruff>=0.6",
    "black>=24.0",
    "sphinx>=5.0.0",
    "myst-parser>=2.0",  # Markdown docs in the Sphinx toctree
    "mypy>=1.0.0",
    "types-PyYAML>=6.0.12",
    "build>=1.2",
    "twine>=5.0",
    "bump2version>=1.0",
]

setup(
    author="Jack Tierney",
    author_email="jackcurragh@gmail.com",
    # Align supported versions with CI and container images
    python_requires=">=3.10,<3.13",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    description="A python command-line utility for the generation of comprehensive reports on the quality of ribosome profiling (Ribo-Seq) datasets",
    entry_points={
        "console_scripts": [
            "RiboMetric=RiboMetric.cli:main",
            "ribometric=RiboMetric.cli:main",
        ],
    },
    install_requires=requirements,
    extras_require={
        "pdf": pdf_requirements,
        "test": test_requirements,
        "dev": dev_requirements,
    },
    license="MIT license",
    long_description=readme,
    long_description_content_type="text/x-rst",
    include_package_data=True,
    keywords="RiboMetric",
    # Lowercase deliberately (REPO_CONTRACT §5.3). PyPI rejects a wheel whose
    # filename is not the normalised project name, with a 400 at upload that
    # `twine check` does not catch -- three TranslonScorer releases died on
    # exactly this. Declaring the name lowercase makes the built filenames
    # `ribometric-*` regardless of the setuptools version; the CI build job
    # asserts it as well. The import package stays `RiboMetric`, and both
    # console scripts are unchanged.
    name="ribometric",
    packages=find_packages(
        include=["RiboMetric", "RiboMetric.*"],
        exclude=[
            "sample_data/*",
        ],
    ),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/JackCurragh/RiboMetric",
    version=version,
    zip_safe=False,
)
