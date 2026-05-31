=======
Results
=======

Output formats
--------------

RiboMetric produces several output formats depending on the flags passed:

**HTML** (default)
    Interactive report with plots for each QC metric.  Open in any browser.

**JSON** (default)
    Machine-readable file containing every metric value, the full metagene
    profiles, and the config used.  Use ``--json`` to request JSON without
    HTML, or omit flags to get both.

    Key top-level sections:

    - ``metrics`` — scalar QC scores (periodicity, coverage, alignment stats, etc.)
    - ``alignment_stats`` — duplicate rate, multimapper rate, soft-clip rate,
      total/mapped/unmapped read counts
    - ``read_length_distribution`` — read count per length
    - ``metagene_profile`` — per-position density around start/stop codons
    - ``rust`` — RUST codon metagene and KL divergence (only when FASTA provided
      and ``rust_mean_kl_divergence`` enabled)
    - ``config`` — the full configuration used for this run

**Summary TSV** (``--summary-tsv``)
    One row per sample, one column per metric.  Convenient for multi-sample
    comparisons in a spreadsheet or downstream script.

**QC Status JSON** (``--qc-status``)
    Machine-readable pass/warn/fail for each metric, using the built-in or
    user-supplied thresholds.

**Comparison CSV** (``--comparison-csv``)
    Wide-format CSV suitable for multi-sample comparison tables.

**Metrics Table CSV** (``--metrics-table``)
    Detailed metrics with per-read-length breakdowns.

**PDF**
    Static version of the HTML report for archiving.  Requires the
    ``RiboMetric[pdf]`` extras: ``pip install RiboMetric[pdf]``.

Pipeline outputs
----------------

Use ``--improved-outputs`` to write all pipeline-friendly formats in one go::

    RiboMetric run -b sample.bam -a annotation.tsv --improved-outputs

Or combine specific flags::

    RiboMetric run -b sample.bam -a annotation.tsv \
        --summary-tsv \
        --qc-status \
        --metrics-table

Example output files can be found in the
`example-reports <https://github.com/JackCurragh/RiboMetric/tree/main/example-reports>`_
directory on GitHub.

See also
--------

:doc:`usage` for full CLI reference.

``REPORTING_GUIDE.md`` in the repository root for format details and examples.
