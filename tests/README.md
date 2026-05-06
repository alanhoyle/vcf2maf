# Tests

This directory is recreated from the upstream `../vcf2maf/tests` fixtures and
golden outputs. The executable test harness is
`tests/test_upstream_docker.py`, a Python `unittest` port of the upstream
Perl `.t` files.

Build the local Docker image first:

```bash
docker build -t vcf2maf:main .
```

## Quick Reference Setup

Use the helper script to download and unpack the small reference data needed by
the Docker-backed tests:

```bash
tests/download_references.sh
```

By default, the script installs:

- `tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa`
- `tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa.fai`
- the VEP 112 GRCh38 chr21 cache under `tests/homo_sapiens/`

To also install the large full GRCh37 VEP cache required by the upstream
`maf2maf` tests, run:

```bash
tests/download_references.sh --full-grch37
```

The script is idempotent and skips files/directories that already exist. Set
`VCF2MAF_DOCKER_IMAGE` or `VCF2MAF_DOCKER_PLATFORM` if you use a different
Docker image tag or platform.

## Reference FASTA

The VEP-backed tests need the chr21 GRCh38 FASTA in this directory:

```bash
curl -LO https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.21.fa tests/
docker run --rm -v "$PWD/tests:/opt/tests" vcf2maf:main samtools faidx tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa
```

## VEP Cache

The tests run VEP with `--vep-data tests`, so the cache must extract to
`tests/homo_sapiens/...`. The upstream comparison suite uses a small chr21-only
VEP 112 GRCh38 cache:

```bash
curl -L -o tests/homo_sapiens_vep_112_GRCh38_chr21.tar.gz https://data.cyri.ac/homo_sapiens_vep_112_GRCh38_chr21.tar.gz
tar -zxf tests/homo_sapiens_vep_112_GRCh38_chr21.tar.gz -C tests/
```

After extraction, this should exist:

```bash
ls tests/homo_sapiens
```

The cache archive and extracted `tests/homo_sapiens/` directory are ignored by
git.

The chr21 cache is enough for the `vcf2maf` tests that use
`tests/test_b38.vcf`. The upstream `maf2maf` fixtures include GRCh37 variants
outside chr21, so those tests skip unless a compatible full GRCh37 VEP cache is
available under `tests/homo_sapiens/112_GRCh37`.

## Full GRCh37 VEP Cache

To run the upstream `maf2maf` tests, install the full homo sapiens VEP 112
GRCh37 cache into `tests/`. The local Docker image includes `vep_install`, so
you can let VEP download and unpack the cache in the expected layout:

```bash
docker run --rm -it \
  --platform linux/amd64 \
  -v "$PWD/tests:/opt/tests" \
  vcf2maf:main \
  vep_install \
    -a cf \
    -s homo_sapiens \
    -y GRCh37 \
    -c /opt/tests \
    --CACHE_VERSION 112
```

After installation, this directory should exist:

```bash
ls tests/homo_sapiens/112_GRCh37
```

This cache is large and is ignored by git via `tests/homo_sapiens/`.

Run the suite:

```bash
python3 -m unittest discover -s tests -v
```

Set `VCF2MAF_DOCKER_IMAGE` to test a differently tagged image.
Set `VCF2MAF_DOCKER_PLATFORM` if you build for something other than
`linux/amd64`.
