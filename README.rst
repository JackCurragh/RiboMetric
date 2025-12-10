================
RiboMetric
================


.. image:: https://img.shields.io/pypi/v/RiboMetric.svg
        :target: https://pypi.python.org/pypi/RiboMetric
        :alt: PyPI version

.. image:: https://readthedocs.org/projects/RiboMetric/badge/?version=latest
        :target: https://RiboMetric.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

.. image:: https://github.com/JackCurragh/RiboMetric/workflows/Build%20and%20Deploy%20Package/badge.svg
        :target: https://github.com/JackCurragh/RiboMetric/actions
        :alt: Build Status

.. image:: https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10-blue
        :target: https://www.python.org/downloads/
        :alt: Python Version

.. image:: https://img.shields.io/github/license/JackCurragh/RiboMetric
        :target: https://github.com/JackCurragh/RiboMetric/blob/main/LICENSE
        :alt: License


A python command-line utility for the generation of comprehensive reports on the quality of ribosome profiling (Ribo-Seq) datasets 


* Free software: MIT license
* Documentation: https://RiboMetric.readthedocs.io.

Installation
------------

To install RiboMetric:

.. code-block:: console

    $ pip install RiboMetric

For PDF export support (adds ~30 dependencies):

.. code-block:: console

    $ pip install RiboMetric[pdf]

Usage
------------

Create annotation files from gff files:

.. code-block:: console

    $ RiboMetric prepare -g gff_file.gff

Use the annotation file to run RiboMetric on a bam file:

.. code-block:: console

    $ RiboMetric run -b bam_file.bam -a annotation_RiboMetric.tsv

View results interactively in your terminal:

.. code-block:: console

    $ RiboMetric view output_RiboMetric_data.json

By default, RiboMetric calculates standard Ribo-Seq QC metrics. To enable optional (theoretical) metrics:

.. code-block:: console

    $ RiboMetric run -b bam_file.bam -a annotation_RiboMetric.tsv --enable-optional-metrics

Or enable specific metrics:

.. code-block:: console

    $ RiboMetric run -b bam_file.bam -a annotation_RiboMetric.tsv --enable-metric periodicity_fourier

For more information on how to use RiboMetric, see the documentation_ or use :code:`--help`

.. _documentation: https://ribometric.readthedocs.io/en/latest/?version=latest

Features
--------

RiboMetric calculates comprehensive quality metrics for Ribo-Seq data:

**Default Metrics (Standard Ribo-Seq QC):**
  * Read length distribution (IQR, coefficient of variation, max proportion)
  * Terminal nucleotide bias (5' and 3' ligation bias detection)
  * 3-nt periodicity (frame dominance and information content)
  * Metagene uniformity (entropy-based)
  * CDS coverage
  * Regional distribution (5'UTR, CDS, 3'UTR proportions and ratios)

**Optional Metrics (Theoretical/Experimental):**
  * Alternative periodicity methods (autocorrelation, Fourier transform, Trips-Viz)
  * Alternative uniformity methods (autocorrelation, Theil index, Gini index)
  * Additional read length metrics (bimodality, normality tests)

Use ``--enable-optional-metrics`` to calculate all metrics, or ``--enable-metric <name>`` for specific ones.

Output Formats
--------------

RiboMetric provides multiple output formats for different use cases:

**For Pipeline Integration:**
  * Summary TSV - One-line summary per sample for quick QC decisions
  * QC Status JSON - Machine-readable pass/warn/fail with thresholds
  * Comparison CSV - Wide format for multi-sample comparison

**For Sample Review:**
  * Interactive TUI - Terminal-based viewer for exploring metrics (``RiboMetric view``)
  * Interactive HTML - Professional reports with executive summary and searchable metrics
  * PDF - Archivable reports for documentation
  * Metrics Table CSV - Detailed metrics with read-length breakdowns

See REPORTING_GUIDE.md_ for complete documentation and examples.

.. _REPORTING_GUIDE.md: REPORTING_GUIDE.md

Requirements
------------

  * Transcriptomic alignments are required in BAM format
  * GFF annotations from Ensembl are also required

Testing
-------

RiboMetric has a comprehensive test suite to ensure reliability:

.. code-block:: console

    $ pip install -r requirements_test.txt
    $ pytest

For more information, see TESTING.md_

.. _TESTING.md: TESTING.md

Credits
-------

This project was worked on by `Lukas Wierdsma`_ during his `Internship at the UCC`_ for Bioinformatics, Howest in 2023.

.. _`Lukas Wierdsma`: https://github.com/Lukas-Wierdsma
.. _`Internship at the UCC`: https://github.com/Lukas-Wierdsma/Internship-UCC-2023/wiki

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
