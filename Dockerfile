FROM debian:bookworm-slim AS builder
ARG TARGETARCH

# Install debootstrap and create a minimal base install in /install_root
RUN apt-get update && \
    apt-get install -y debootstrap ca-certificates curl && \
    debootstrap --variant=minbase --include=ca-certificates bookworm /install_root

# Download and install conda into the builder stage's /opt/conda
ENV MINICONDA_VERSION=py312_25.5.1-1
RUN set -eux; \
    arch="${TARGETARCH:-$(uname -m)}"; \
    case "$arch" in \
        amd64|x86_64) miniconda_arch="x86_64" ;; \
        arm64|aarch64) miniconda_arch="aarch64" ;; \
        *) echo "Unsupported Docker target architecture: $arch" >&2; exit 1 ;; \
    esac; \
    curl -sL "https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-${miniconda_arch}.sh" -o /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -bup /opt/conda && \
    rm -f /tmp/miniconda.sh

# Modify PATH to see conda, use faster dependency solver, and accept TOS
ENV PATH="/opt/conda/bin:${PATH}"
RUN conda config --set solver libmamba && \
    conda tos accept --channel defaults

# Use conda to install vcf2maf tools/dependencies into /usr/local
ENV VEP_VERSION=112.0 \
    HTSLIB_VERSION=1.20 \
    BCFTOOLS_VERSION=1.20 \
    SAMTOOLS_VERSION=1.20 \
    LIFTOVER_VERSION=447
RUN conda create -y -p /usr/local && \
    conda install -y -p /usr/local \
    -c conda-forge \
    -c bioconda \
    -c defaults \
    ensembl-vep==${VEP_VERSION} \
    perl-list-moreutils \
    htslib==${HTSLIB_VERSION} \
    bcftools==${BCFTOOLS_VERSION} \
    samtools==${SAMTOOLS_VERSION} \
    ucsc-liftover==${LIFTOVER_VERSION}

# Install the vcf2maf Python package into the conda environment at /usr/local
COPY pyproject.toml README.md /tmp/vcf2maf/
COPY vcf2maf/ /tmp/vcf2maf/vcf2maf/
RUN /usr/local/bin/pip install --no-build-isolation /tmp/vcf2maf/

# Deploy the minimal OS and tools into a clean target layer
FROM scratch

LABEL maintainer="Cyriac Kandoth <ckandoth@gmail.com>"

COPY --from=builder /install_root /
COPY --from=builder /usr/local /usr/local
COPY data /opt/data
WORKDIR /opt
