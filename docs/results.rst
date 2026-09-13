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

    The file has two top-level keys, ``results`` and ``config`` (the full
    configuration used for the run). Inside ``results``:

    - ``metrics`` — raw measurements, in natural units and natural direction
    - ``scores`` — 0–1 scores, always higher-is-better, named for the good property
    - ``metrics_legacy`` — every pre-2.0 key with its pre-2.0 value; removed in 2.1
    - ``provenance`` — input file hashes and the effective configuration hash
    - ``offsets`` — the A-site offsets applied, per read length
    - ``alignment_stats`` — duplicate, multimapper and soft-clip rates; read counts
    - ``read_length_distribution`` — read count per length
    - ``metagene_profile`` — per-position density around start and stop codons
    - ``rust`` — RUST codon metagene and KL divergence (only with a FASTA and
      ``--enable-metric rust_mean_kl_divergence``)

    The split between ``metrics`` and ``scores`` is new in 2.0; see
    :doc:`METRIC_NAMING`.

**Summary TSV** (``--summary-tsv``)
    One row per sample, one column per metric.  Convenient for multi-sample
    comparisons in a spreadsheet or downstream script.

**QC Status JSON** (written by default)
    Machine-readable pass/warn/fail verdict. It comes from the same scores and
    thresholds as the HTML report, including any ``scoring:`` overrides in your
    config, so the two always agree.

**Comparison CSV** (``--comparison-csv``)
    Wide-format CSV suitable for multi-sample comparison tables.

**Metrics Table CSV** (``--metrics-table``)
    Detailed metrics with per-read-length breakdowns.

**PDF**
    Static version of the HTML report for archiving.  Requires the
    ``pdf`` extra: ``pip install "ribometric[pdf]"``.

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

:doc:`REPORTING_GUIDE` for format details and worked examples.
