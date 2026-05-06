#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TESTS_DIR="${ROOT_DIR}/tests"

IMAGE="${VCF2MAF_DOCKER_IMAGE:-vcf2maf:main}"
PLATFORM="${VCF2MAF_DOCKER_PLATFORM:-linux/amd64}"

GRCH38_FASTA="Homo_sapiens.GRCh38.dna.chromosome.21.fa"
GRCH38_FASTA_GZ="${GRCH38_FASTA}.gz"
GRCH38_FASTA_URL="https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/${GRCH38_FASTA_GZ}"

GRCH38_CACHE_ARCHIVE="homo_sapiens_vep_112_GRCh38_chr21.tar.gz"
GRCH38_CACHE_URL="https://data.cyri.ac/${GRCH38_CACHE_ARCHIVE}"

FULL_GRCH37=0

usage() {
  cat <<EOF
Usage: $0 [--full-grch37]

Downloads test reference data into ${TESTS_DIR}.

Default:
  - GRCh38 chromosome 21 FASTA
  - samtools FASTA index
  - VEP 112 GRCh38 chr21 cache

Optional:
  --full-grch37  Also install the large full VEP 112 GRCh37 cache for maf2maf tests

Environment:
  VCF2MAF_DOCKER_IMAGE     Docker image to use for samtools/vep_install [${IMAGE}]
  VCF2MAF_DOCKER_PLATFORM  Docker platform [${PLATFORM}]
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --full-grch37)
      FULL_GRCH37=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

require_command() {
  local cmd="$1"
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "ERROR: Required command not found: ${cmd}" >&2
    exit 1
  fi
}

require_docker_image() {
  if ! docker image inspect "${IMAGE}" >/dev/null 2>&1; then
    echo "ERROR: Docker image '${IMAGE}' not found." >&2
    echo "Build it first, for example: docker build -t ${IMAGE} ." >&2
    exit 1
  fi
}

download_file() {
  local url="$1"
  local output="$2"
  if [[ -s "${output}" ]]; then
    echo "Already present: ${output}"
    return
  fi
  echo "Downloading ${url}"
  curl -L --fail --retry 3 --output "${output}" "${url}"
}

index_fasta() {
  local fasta="${TESTS_DIR}/${GRCH38_FASTA}"
  if [[ -s "${fasta}.fai" ]]; then
    echo "Already indexed: ${fasta}.fai"
    return
  fi
  require_docker_image
  echo "Indexing ${fasta}"
  docker run --rm \
    --platform "${PLATFORM}" \
    -v "${TESTS_DIR}:/opt/tests" \
    "${IMAGE}" \
    samtools faidx "/opt/tests/${GRCH38_FASTA}"
}

install_grch38_fasta() {
  local fasta="${TESTS_DIR}/${GRCH38_FASTA}"
  local fasta_gz="${TESTS_DIR}/${GRCH38_FASTA_GZ}"

  if [[ ! -s "${fasta}" ]]; then
    download_file "${GRCH38_FASTA_URL}" "${fasta_gz}"
    echo "Decompressing ${fasta_gz}"
    gzip -dkf "${fasta_gz}"
  else
    echo "Already present: ${fasta}"
  fi

  index_fasta
}

install_grch38_cache() {
  local archive="${TESTS_DIR}/${GRCH38_CACHE_ARCHIVE}"
  if [[ -d "${TESTS_DIR}/homo_sapiens/112_GRCh38" ]]; then
    echo "Already present: ${TESTS_DIR}/homo_sapiens/112_GRCh38"
    return
  fi
  download_file "${GRCH38_CACHE_URL}" "${archive}"
  echo "Extracting ${archive}"
  tar -zxf "${archive}" -C "${TESTS_DIR}"
}

install_full_grch37_cache() {
  if [[ -d "${TESTS_DIR}/homo_sapiens/112_GRCh37" ]]; then
    echo "Already present: ${TESTS_DIR}/homo_sapiens/112_GRCh37"
    return
  fi
  require_docker_image
  echo "Installing full VEP 112 GRCh37 cache. This is large and can take a while."
  docker run --rm -it \
    --platform "${PLATFORM}" \
    -v "${TESTS_DIR}:/opt/tests" \
    "${IMAGE}" \
    vep_install \
      -a cf \
      -s homo_sapiens \
      -y GRCh37 \
      -c /opt/tests \
      --CACHE_VERSION 112
}

require_command curl
require_command gzip
require_command tar
require_command docker

mkdir -p "${TESTS_DIR}"

install_grch38_fasta
install_grch38_cache

if [[ "${FULL_GRCH37}" -eq 1 ]]; then
  install_full_grch37_cache
else
  cat <<EOF

Skipping full GRCh37 VEP cache.
Run this to enable the maf2maf upstream tests:

  $0 --full-grch37

EOF
fi

echo "Reference setup complete."
