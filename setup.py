#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

# Core dependencies needed to run RiboMetric
requirements = [
    "biopython>=1.81",
    "click>=8.0",
    "gffpandas>=1.2.0",
    "Jinja2>=3.0",
    "kaleido==0.2.1",
    "numpy>=1.24",
    "oxbow>=0.3.0",
    "pandas>=2.0",
    "plotly>=5.14",
    "pyarrow>=10.0.0",
    "pysam>=0.21.0",
    "PyYAML>=6.0",
    "rich>=13.0.0",
    "scipy>=1.10",
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
dev_requirements = test_requirements + [
    "flake8>=6.0.0",
    "black>=23.0.0",
    "sphinx>=5.0.0",
]

setup(
    author="Jack Tierney",
    author_email="jackcurragh@gmail.com",
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
    ],
    description="A python command-line utility for the generation of comprehensive reports on the quality of ribosome profiling (Ribo-Seq) datasets",
    entry_points={
        "console_scripts": [
            "RiboMetric=RiboMetric.cli:main",
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
    name="RiboMetric",
    packages=find_packages(
        include=["RiboMetric", "RiboMetric.*"],
        exclude=[
            "sample_data/*",
        ],
    ),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/JackCurragh/RiboMetric",
    version="0.1.12",
    zip_safe=False,
)
