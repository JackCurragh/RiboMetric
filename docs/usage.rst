.. role:: bash(code)
   :language: bash
..
    Lukas: Draft for usage docs, feedback appreciated :)
=====
Usage
=====
..
    Lukas: is this still relevant?
To use RiboMetric in a project::

    import RiboMetric

After installing RiboMetric, you can use the command line application to run RiboMetric on bam files with::

    RiboMetric run -b bam_file.bam
    
This will return a basic report containing metrics for ligation bias, read frame periodicity and more.
By adding an anottation file to the command, RiboMetric will be able to generate more metrics with the extra information.
These annotation files are created from gff files by using the command::

    RiboMetric preprare -g gff_file.gff
    
Once the annotation file is created, run the bam file with the annotation file::

    RiboMetric run -b bam_file.bam -a annotation_RiboMetric.tsv
    
This will add on metrics for metagene profile and mRNA distribution. Further options for the command can be called with :bash:`--help`::

    --subsample: Number of reads to subsample from the bam file (default = 1000000)
    --json: Output the results as a json file
    --...
