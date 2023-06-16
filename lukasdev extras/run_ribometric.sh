#!/bin/bash

# Directory containing the bam files
bam_dir="/home/DATA2/jack/Human_Transcriptome_Alignments/bams"

# Annotation file
annotation_file="/home/lukas/projects/RiboMetric/sample_data/gencode.v25.annotation_RiboMetric.tsv"

# Iterate over each bam file in the directory
for  bam_file in "$bam_dir"/*.bam_sorted; do
	# Copy the bam file to the current directory
	echo "cp $bam_file ."
	cp "$bam_file" .

	# Get the filename of the copied bam file
	temp_bam=$(basename "$bam_file")

	# Run samtools index on copied bam file
	echo "samtools index $temp_bam"
	samtools index "$temp_bam"

	# Run RiboMetric command for copied bam file with memory profiler
	echo "mprof run RiboMetric run -b $temp_bam -a $annotation_file" --all -p 12
	mprof run RiboMetric run -b "$temp_bam" -a "$annotation_file" --all -p 12

	# Export the mprof plots
	mprof plot -o mprof/"$temp_bam"

	# Remove mprof data
	mprof clean

	# Delete the temp files
	echo "rm $temp_bam ${temp_bam}.bai"
	rm "$temp_bam" "${temp_bam}.bai"

	# Print a splitting line
	echo "---------------------------------------------------"
done

