.. highlight:: shell

============
Installation
============

From PyPI (recommended)
------------------------

.. code-block:: console

    $ pip install ribometric

This installs the latest release (the name is case-insensitive, so
``pip install RiboMetric`` works too).  For PDF export support:

.. code-block:: console

    $ pip install "ribometric[pdf]"

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

Images are published to the GitHub Container Registry:

* ``:latest`` — the most recent release.
* ``:vX.Y.Z`` and ``:X.Y`` — a specific release; pin one of these for
  reproducible runs. (Releases up to 1.4.3 were tagged without the ``v``, for
  example ``:1.4.3``.)
* ``:dev`` — the tip of the ``dev`` branch: the newest features, not yet
  released.

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
