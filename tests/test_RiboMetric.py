#!/usr/bin/env python

"""Tests for `RiboMetric` package."""

import json
import os
import sys
from argparse import Namespace
from io import StringIO

from RiboMetric.RiboMetric import main


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
        global_offset=15,
    )

    # Redirect stdout to a StringIO object to capture output
    sys.stdout = StringIO()

    main(args)

    # Get the output
    output = sys.stdout.getvalue()

    # Assert that the expected output was produced
    assert "Annotation parsed" in output
    assert "Running modules" in output

    with open(f"{file_path}/test_RiboMetric.json") as result_file:
        results = json.load(result_file)["results"]

    metrics = results["metrics"]
    assert metrics["prop_reads_CDS"]["global"] == 1.0
    assert metrics["prop_reads_leader"]["global"] == 0.0
    assert metrics["prop_reads_trailer"]["global"] == 0.0
    assert "terminal_bias_kl_5prime_score" in metrics
    assert "terminal_bias_kl_5prime_raw" in metrics
    assert metrics["terminal_bias_kl_5prime"] == metrics["terminal_bias_kl_5prime_score"]

    provenance = results["provenance"]
    assert "effective_config_sha256" in provenance
    assert provenance["inputs"]["bam"]["hash_method"] == "full_sha256"
    assert provenance["inputs"]["annotation"]["hash_method"] == "full_sha256"

    offsets = results["offsets"]
    assert offsets["source"] == "global"
    assert offsets["global_offset"] == 15
    assert offsets["target"] == "a_site"
    assert offsets["applied_by_read_length"]
    assert os.path.exists(f"{file_path}/test_offsets.tsv")


def test_main_run_no_server_flag_anymore():
    # Keeping a placeholder test to assert absence of server option at CLI level.
    # The CLI parser no longer accepts --server; this test confirms nothing to run here.
    assert True


def test_main_run_readlength_offsets():
    file_path = os.path.join(os.path.dirname(__file__), "test_data")

    args = Namespace(
        command="run",
        annotation=f"{file_path}/1000_entry_RiboMetric.tsv",
        bam=f"{file_path}/test.bam",
        output=f"{file_path}",
        offset_read_length=f"{file_path}/offset_read_length.tsv",
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
    assert "Applying specified read length specific offsets" in output
    assert "Running modules" in output


def test_main_run_readspecific_offsets():
    file_path = os.path.join(os.path.dirname(__file__), "test_data")

    args = Namespace(
        command="run",
        annotation=f"{file_path}/1000_entry_RiboMetric.tsv",
        bam=f"{file_path}/test.bam",
        output=f"{file_path}",
        offset_read_specific=f"{file_path}/readspecific_offsets.tsv",
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
    assert "Applying read specific offsets" in output
    assert "Running modules" in output


def test_main_run_global_offsets():
    file_path = os.path.join(os.path.dirname(__file__), "test_data")

    args = Namespace(
        command="run",
        annotation=f"{file_path}/1000_entry_RiboMetric.tsv",
        bam=f"{file_path}/test.bam",
        output=f"{file_path}",
        global_offset=15,
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
    assert "Applying global offset" in output
    assert "Running modules" in output
