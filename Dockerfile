FROM python:3.10-slim

LABEL description="RiboMetric container for ribosome profiling QC" \
      author="Jack Curragh" \
      org.opencontainers.image.source="https://github.com/JackCurragh/RiboMetric"

# System deps required to build pysam and run samtools
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    samtools \
    libbz2-dev \
    zlib1g-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . .

RUN pip install .

WORKDIR /data

CMD ["/bin/bash"]
