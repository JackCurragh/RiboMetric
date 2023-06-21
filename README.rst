================
RiboMetric
================


.. image:: https://img.shields.io/pypi/v/RiboMetric.svg
        :target: https://pypi.python.org/pypi/RiboMetric

.. image:: https://readthedocs.org/projects/RiboMetric/badge/?version=latest
        :target: https://RiboMetric.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status

.. image:: https://pyup.io/repos/github/JackCurragh/RiboMetric/shield.svg
     :target: https://pyup.io/repos/github/JackCurragh/RiboMetric/
     :alt: Updates


A python command-line utility for the generation of comprehensive reports on the quality of ribosome profiling (Ribo-Seq) datasets 


* Free software: MIT license
* Documentation: https://RiboMetric.readthedocs.io.

Installation
------------

To install RiboMetric:

.. code-block:: console

    $ pip install RiboMetric

Usage
------------

Create annotation files from gff files:

.. code-block:: console

    $ RiboMetric preprare -g gff_file.gff
    
Use the annotation file to run RiboMetric on a bam file:

.. code-block:: console

    $ RiboMetric run -b bam_file.bam -a annotation_RiboMetric.tsv

For more information on how to use RiboMetric, see the documentation_ or use :code:`--help`

.. _documentation: https://ribometric.readthedocs.io/en/latest/?version=latest

Features
--------

There are a number of limitations to running RiboMetric at the moment. These include:

  * Transciptomic alignments are required in BAM format 
  * GFF annotations from Ensembl are also requried

Credits
-------

This project was worked on by `Lukas Wierdsma`_ during his `Internship at the UCC`_ for Bioinformatics, Howest in 2023.

.. _`Lukas Wierdsma`: https://github.com/Lukas-Wierdsma
.. _`Internship at the UCC`: https://github.com/Lukas-Wierdsma/Internship-UCC-2023/wiki

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
