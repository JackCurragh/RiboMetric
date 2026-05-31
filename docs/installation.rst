.. highlight:: shell

============
Installation
============

From PyPI (recommended)
------------------------

.. code-block:: console

    $ pip install RiboMetric

This installs the latest stable release.  For PDF export support:

.. code-block:: console

    $ pip install RiboMetric[pdf]

From source
-----------

Clone the repository and install in editable mode:

.. code-block:: console

    $ git clone https://github.com/JackCurragh/RiboMetric
    $ cd RiboMetric
    $ pip install -e .

For development (includes test and lint dependencies):

.. code-block:: console

    $ pip install -e ".[dev]"

Using pixi
----------

`pixi <https://pixi.sh>`_ resolves all dependencies including test extras:

.. code-block:: console

    $ git clone https://github.com/JackCurragh/RiboMetric
    $ cd RiboMetric
    $ pixi run test   # installs everything and runs the test suite

Docker
------

A Docker image is built and pushed on every tagged release:

.. code-block:: console

    $ docker pull ghcr.io/jackcurragh/ribometric:latest
    $ docker run --rm -v $(pwd):/data \
        ghcr.io/jackcurragh/ribometric:latest \
        RiboMetric run -b /data/sample.bam

Requirements
------------

* Python 3.10, 3.11, or 3.12
* ``samtools`` >= 1.10 on your PATH (used for BAM indexing and idxstats)
* A coordinate-sorted, indexed BAM file aligned to a transcriptome reference
