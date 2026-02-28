FROM python:3.10

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

# Avoid site‑packages leakage and disable static image export by default.
# Plotly's Kaleido requires Chrome; skipping images avoids runtime errors in
# headless containers unless the user explicitly opts in to PDF/static PNGs.
ENV PYTHONNOUSERSITE=1 \
    RIBOMETRIC_SKIP_IMAGES=1

WORKDIR /app

COPY . .

# Install package from the repo snapshot copied into the image.
# If you need PDF/static images, install extras and Chrome at build time:
#   RUN pip install -U plotly kaleido 'RiboMetric[pdf]' && \
#       python -c 'import kaleido; kaleido.get_chrome_sync()' && \
#       echo "RIBOMETRIC_SKIP_IMAGES=0" >> /etc/environment
RUN pip install .

WORKDIR /data

CMD ["/bin/bash"]
