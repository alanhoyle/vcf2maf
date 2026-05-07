# Tests

This directory contains the upstream fixture inputs, golden MAF outputs from
the Perl implementation, and the Python `unittest` test harness
`test_upstream_docker.py`.

The suite runs the four tools inside Docker and diffs the output against the
upstream Perl golden MAFs.

---

## 1. Build the Docker image

```bash
docker build -t vcf2maf:main .
```

---

## 2. Download reference data

Use the helper script to fetch and index the reference data needed by the tests:

```bash
tests/download_references.sh
```

By default this installs:

- `tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa` (+ `.fai` index)
- VEP 112 GRCh38 chr21 cache under `tests/homo_sapiens/`

The script is idempotent and skips files that already exist.

### Full GRCh37 cache (required for `maf2maf` tests)

The upstream `maf2maf` fixtures include GRCh37 variants outside chr21, so those
tests are automatically skipped unless a full GRCh37 VEP 112 cache is present
under `tests/homo_sapiens/112_GRCh37`.  To install it (~20 GB):

```bash
tests/download_references.sh --full-grch37
```

Or install manually via the VEP installer bundled in the Docker image:

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

After installation this directory should exist:

```bash
ls tests/homo_sapiens/112_GRCh37
```

The cache directories are ignored by git.

---

## 3. Run the suite

```bash
python3 -m unittest discover -s tests -v
```

The `maf2maf` tests are automatically skipped when the full GRCh37 cache is absent.

---

## Environment variables

| Variable                  | Purpose                                                      |
| ------------------------- | ------------------------------------------------------------ |
| `VCF2MAF_DOCKER_IMAGE`    | Docker image to test against (default: `vcf2maf:main`)       |
| `VCF2MAF_DOCKER_PLATFORM` | Docker platform (default: `linux/amd64`)                     |
| `PRESERVE_TESTS`          | Set to `1` to keep intermediate output files after test runs |

---

## Reference data layout

```text
tests/
├── homo_sapiens/              # VEP cache (git-ignored)
│   ├── 112_GRCh38/            # chr21-only cache (small, required)
│   └── 112_GRCh37/            # full cache (large, optional)
├── Homo_sapiens.GRCh38.dna.chromosome.21.fa      # reference FASTA
├── Homo_sapiens.GRCh38.dna.chromosome.21.fa.fai  # samtools index
├── test.vcf / test.maf / test_b38.vcf / …        # upstream fixture inputs
└── *.maf                                          # Perl golden outputs
```

### Manual FASTA setup

If you need to fetch and index the reference FASTA by hand:

```bash
curl -LO https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
gunzip Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
mv Homo_sapiens.GRCh38.dna.chromosome.21.fa tests/
docker run --rm -v "$PWD/tests:/opt/tests" vcf2maf:main \
  samtools faidx /opt/tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa
```
