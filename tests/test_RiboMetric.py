#!/usr/bin/env python

"""Tests for `RiboMetric` package."""

from RiboMetric.RiboMetric import main

from argparse import Namespace
from io import StringIO
import sys
import os


def test_main_prepare():
    file_path = os.path.join(os.path.dirname(__file__), "test_data")

    args = Namespace(
        command="prepare",
        gff=f"{file_path}/1000_entry.gff",
        transcripts=1000,
        output=f"{file_path}",
        config=f"{file_path}/../../config.yml",
        threads=2,
    )

    # Redirect stdout to a StringIO object to capture output
    sys.stdout = StringIO()

    main(args)

    # Get the output
    output = sys.stdout.getvalue()

    # Assert that the expected output was produced
    assert "Parsing gff" in output


def test_main_run():
    file_path = os.path.join(os.path.dirname(__file__), "test_data")

    args = Namespace(
        command="run",
        annotation=f"{file_path}/1000_entry_RiboMetric.tsv",
        bam=f"{file_path}/test.bam",
        output=f"{file_path}",
        config=f"{file_path}/../../config.yml",
        all=False,
        gff=None,
        fasta=None,
        subsample=1000,
        transcripts=None,
        json=True,
        html=True,
        pdf=None,
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


def test_main_run_server():
    file_path = os.path.join(os.path.dirname(__file__), "test_data")
    args = Namespace(
        command="run",
        annotation=f"{file_path}/1000_entry_RiboMetric.tsv",
        bam=f"{file_path}/test.bam",
        output=f"{file_path}",
        config=f"{file_path}/../../config.yml",
        all=False,
        gff=None,
        fasta=None,
        subsample=1000,
        transcripts=None,
        json=None,
        html=None,
        pdf=None,
        csv=True,
        server=True
    )

    # Redirect stdout to a StringIO object to capture output
    sys.stdout = StringIO()

    main(args)

    # Get the output
    output = sys.stdout.getvalue()

    # Assert that the expected output was produced
    assert 'Annotation parsed' in output
    assert 'Running modules' in output


def test_main_run_json():
    file_path = os.path.join(os.path.dirname(__file__), "test_data")
    print(file_path)
    args = Namespace(
        command="run",
        json_in=f"{file_path}/test.json",
        output=f"{file_path}",
        config=f"{file_path}/../../config.yml",
        all=False,
        gff=None,
        fasta=None,
        subsample=1000,
        transcripts=None,
        json=None,
        html=None,
        pdf=True,
        csv=None,
    )

    # Redirect stdout to a StringIO object to capture output
    sys.stdout = StringIO()

    main(args)

    # Get the output
    output = sys.stdout.getvalue()

    # Assert that the expected output was produced
    assert 'JSON' in output
    assert 'Your pdf report can be found' in output
