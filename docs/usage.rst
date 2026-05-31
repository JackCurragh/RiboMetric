.. highlight:: shell

=====
Usage
=====

Quick start
-----------

Run RiboMetric on a transcriptome-mapped BAM file::

    RiboMetric run -b sample.bam

This produces an HTML report and JSON output with standard read-level QC metrics
(read length distribution, terminal nucleotide bias, periodicity, alignment stats).

Annotation-aware run
--------------------

Create an annotation file from a GFF3/GTF once per genome release::

    RiboMetric prepare -g annotation.gff3 -o outdir/

Then run with the annotation to unlock metagene, CDS coverage, regional
distribution, and codon-level metrics::

    RiboMetric run -b sample.bam -a outdir/annotation_RiboMetric.tsv

Or pass a pre-built annotation TSV directly if you already have one::

    RiboMetric run -b sample.bam -a homo_sapiens_grch38_114_ribometric.tsv

RUST metric (requires a transcriptome FASTA)
--------------------------------------------

The RUST (Ribo-seq Unit Step Transformation) metric quantifies codon-specific
A-site accumulation. It requires a transcriptome FASTA whose transcript
coordinates match the annotation (use the same Ensembl/GENCODE release).
Gzip-compressed FASTA is supported::

    RiboMetric run -b sample.bam \
        -a annotation.tsv \
        -f Homo_sapiens.GRCh38.cdna.all.fa.gz \
        --enable-metric rust_mean_kl_divergence

Output options
--------------

By default RiboMetric writes an HTML report and a JSON file.  Additional
output formats are available::

    # JSON only (skip HTML rendering):
    RiboMetric run -b sample.bam -a annotation.tsv --json

    # Pipeline-friendly outputs (summary TSV, QC-status JSON, metrics CSV):
    RiboMetric run -b sample.bam -a annotation.tsv --improved-outputs

    # Or individually:
    RiboMetric run -b sample.bam -a annotation.tsv \
        --summary-tsv \
        --qc-status \
        --comparison-csv \
        --metrics-table

Optional metrics
----------------

All metrics in the ``optional`` list in ``config.yml`` are skipped unless
you opt in::

    # Enable all optional metrics:
    RiboMetric run -b sample.bam -a annotation.tsv --enable-optional-metrics

    # Enable one specific metric:
    RiboMetric run -b sample.bam -a annotation.tsv \
        --enable-metric periodicity_fourier

See ``docs/METRICS.md`` for a full description of every metric and which list
it belongs to.

Evaluate (pipeline gating)
---------------------------

Gate a pipeline step on QC thresholds — exits 0 (PASS), 1 (WARN), or 2 (FAIL)::

    RiboMetric evaluate -i sample_RiboMetric_data.json

Provide a YAML of custom thresholds (built-in defaults used if omitted)::

    RiboMetric evaluate -i sample_RiboMetric_data.json -e thresholds.yml -o eval.json

Threshold YAML format::

    thresholds:
      periodicity_dominance:        {pass: 0.70, warn: 0.50}
      prop_reads_CDS:               {pass: 0.70, warn: 0.50}
      duplicate_rate:               {pass: 0.20, warn: 0.40}
      stop_codon_readthrough_ratio: {pass: 0.10, warn: 0.20}

Interactive viewer
------------------

Explore a JSON result in an interactive terminal UI::

    RiboMetric view sample_RiboMetric_data.json

Multimapper handling
--------------------

For STAR transcriptome BAMs, frame-sensitive calculations default to
uniquely-mapped reads only (MAPQ = 255)::

    # Default: uniquely-mapped only
    RiboMetric run -b star.bam -a annotation.tsv --multimap-filter unique_only

    # Pre-1.1 behaviour: all primary reads
    RiboMetric run -b star.bam -a annotation.tsv --multimap-filter none

Sequence-less BAMs
------------------

Some BAMs omit the SEQ field.  Use ``--skip-sequence-metrics`` to avoid errors::

    RiboMetric run -b no_seq.bam -a annotation.tsv --skip-sequence-metrics

Number of reads
---------------

Subsample the BAM to speed up exploratory runs (default: 1 000 000)::

    RiboMetric run -b sample.bam --subsample 500000

Full option reference
---------------------

::

    RiboMetric run --help
    RiboMetric prepare --help
    RiboMetric evaluate --help
    RiboMetric view --help
