# Tests

This directory is recreated from the upstream `../vcf2maf/tests` fixtures and
golden outputs. The executable test harness is
`tests/test_upstream_docker.py`, a Python `unittest` port of the upstream
Perl `.t` files.

Build the local Docker image first:

```bash
docker build -t vcf2maf:main .
```

The VEP-backed tests need the chr21 GRCh38 FASTA in this directory:

```bash
curl -LO https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.21.fa tests/
docker run --rm -v "$PWD/tests:/opt/tests" vcf2maf:main samtools faidx tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa
```

Run the suite:

```bash
python3 -m unittest discover -s tests -v
```

Set `VCF2MAF_DOCKER_IMAGE` to test a differently tagged image.
Set `VCF2MAF_DOCKER_PLATFORM` if you build for something other than
`linux/amd64`.
