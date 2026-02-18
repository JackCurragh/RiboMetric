FROM condaforge/mambaforge:latest

LABEL description="RiboMetric container for ribosome profiling QC" \
      author="Jack Curragh" \
      org.opencontainers.image.source="https://github.com/JackCurragh/RiboMetric"

# Install conda packages (precompiled, avoids compilation issues)
RUN mamba install -y -c conda-forge -c bioconda \
    python=3.10 \
    pip \
    biopython \
    pysam \
    rich \
    pyarrow \
    samtools \
 && mamba clean -a -y

WORKDIR /app

COPY . .

RUN pip install --no-build-isolation .

WORKDIR /data

CMD ["/bin/bash"]
