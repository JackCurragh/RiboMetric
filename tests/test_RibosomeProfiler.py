#!/usr/bin/env python

"""Tests for `RibosomeProfiler` package."""

from RibosomeProfiler.RibosomeProfiler import main

from argparse import Namespace
from io import StringIO
import sys
import os


def test_main_prepare():
    file_path = os.path.join(os.path.dirname(__file__), 'test_data')

    args = Namespace(command='prepare',
                     gff=f'{file_path}/1000_entry.gff',
                     transcripts=1000,
                     output=f'{file_path}',
                     config=f'{file_path}/test_config.txt')

    # Redirect stdout to a StringIO object to capture output
    sys.stdout = StringIO()

    main(args)

    # Get the output
    output = sys.stdout.getvalue()

    # Assert that the expected output was produced
    assert 'Parsing gff' in output


def test_main_run():
    file_path = os.path.join(os.path.dirname(__file__), 'test_data')

    args = Namespace(command='run',
                     annotation=f'{file_path}/1000_entry_RibosomeProfiler.tsv',
                     bam=f'{file_path}/test.bam',
                     output=f'{file_path}',
                     config=f'{file_path}/test_config.txt',
                     all=False,
                     gff=None,
                     fasta=None,
                     subsample=1000,
                     transcripts=None,
                     json=None,
                     html=True,
                     csv=None,
                     )

    # Redirect stdout to a StringIO object to capture output
    sys.stdout = StringIO()

    main(args)

    # Get the output
    output = sys.stdout.getvalue()

    # Assert that the expected output was produced
    assert 'Annotation parsed' in output
    assert 'Running modules' in output


