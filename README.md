# vcf2maf (Python port)

A pure-Python port of [mskcc/vcf2maf](https://github.com/mskcc/vcf2maf), which converts
VCF files into MAF format by annotating each variant to exactly one gene isoform via
[Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html).

---

## Scripts

| Script       | Purpose                                                                    |
| ------------ | -------------------------------------------------------------------------- |
| `vcf2maf.py` | Convert a VCF to MAF (calls VEP, then maps effects to MAF format)          |
| `maf2vcf.py` | Convert a MAF back to per-tumor/normal-pair VCFs                           |
| `maf2maf.py` | Re-annotate an existing MAF (runs `maf2vcf` + `vcf2maf` internally)        |
| `vcf2vcf.py` | Normalise / clean a VCF (multiallelic splitting, left-alignment, liftOver) |

---

## Requirements

- Python ≥ 3.8 (standard library only – no third-party dependencies)
- [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html)
  with the relevant offline cache installed
- `samtools` on `PATH` (used by `maf2vcf` for indel anchor bases)
- `liftOver` on `PATH` (optional, only needed with `--remap-chain`)

---

## Quick start

### VCF → MAF

```bash
python vcf2maf.py \
  --input-vcf  tests/test.vcf \
  --output-maf tests/test.vep.maf \
  --tumor-id   WD1309 \
  --normal-id  NB1308
```

Override VEP/cache paths if not in default locations:

```bash
python vcf2maf.py \
  --input-vcf  tests/test.vcf \
  --output-maf tests/test.vep.maf \
  --tumor-id   WD1309 \
  --normal-id  NB1308 \
  --vep-path   /opt/vep \
  --vep-data   /srv/vep
```

VarScan-style VCF (hardcoded TUMOR/NORMAL column IDs):

```bash
python vcf2maf.py \
  --input-vcf     tests/test_varscan.vcf \
  --output-maf    tests/test_varscan.vep.maf \
  --tumor-id      WD1309 \
  --normal-id     NB1308 \
  --vcf-tumor-id  TUMOR \
  --vcf-normal-id NORMAL
```

Skip VEP (use existing annotations or produce a minimal MAF):

```bash
python vcf2maf.py \
  --input-vcf  tests/test.vep.vcf \
  --output-maf tests/test.maf \
  --inhibit-vep
```

---

### MAF → MAF (re-annotation)

```bash
python maf2maf.py \
  --input-maf  tests/test.maf \
  --output-maf tests/test.vep.maf
```

---

### MAF → VCF

```bash
python maf2vcf.py \
  --input-maf  tests/test.maf \
  --output-dir vcfs/ \
  --ref-fasta  ~/.vep/homo_sapiens/112_GRCh37/Homo_sapiens.GRCh37.dna.toplevel.fa.gz
```

---

### VCF normalisation

```bash
python vcf2vcf.py \
  --input-vcf     input.vcf \
  --output-vcf    normalised.vcf \
  --vcf-tumor-id  TUMOR \
  --vcf-normal-id NORMAL
```

---

## Key options (vcf2maf.py)

| Option            | Default                | Description                                          |
| ----------------- | ---------------------- | ---------------------------------------------------- |
| `--tumor-id`      | `TUMOR`                | Tumor sample barcode in output MAF                   |
| `--normal-id`     | `NORMAL`               | Normal sample barcode in output MAF                  |
| `--vcf-tumor-id`  | `--tumor-id`           | Sample column name in input VCF                      |
| `--vcf-normal-id` | `--normal-id`          | Sample column name in input VCF                      |
| `--vep-path`      | `~/miniconda3/bin`     | Directory containing the `vep` binary                |
| `--vep-data`      | `~/.vep`               | VEP offline cache directory                          |
| `--ref-fasta`     | `~/.vep/…GRCh37…fa.gz` | Reference FASTA (must be samtools-indexed)           |
| `--ncbi-build`    | `GRCh37`               | Genome build; used in MAF `NCBI_Build` column        |
| `--species`       | `homo_sapiens`         | Ensembl species name                                 |
| `--cache-version` | auto                   | VEP offline cache version                            |
| `--inhibit-vep`   | off                    | Skip VEP; parse existing CSQ/ANN if present          |
| `--custom-enst`   | –                      | File of preferred Ensembl transcript IDs             |
| `--retain-info`   | –                      | Comma-sep INFO keys to add as extra MAF columns      |
| `--retain-fmt`    | –                      | Comma-sep FORMAT keys to add as extra MAF columns    |
| `--retain-ann`    | –                      | Comma-sep VEP CSQ fields to add as extra MAF columns |
| `--max-subpop-af` | `0.0004`               | gnomAD AF threshold for `common_variant` FILTER tag  |
| `--min-hom-vaf`   | `0.7`                  | VAF threshold to call a homozygous variant           |
| `--remap-chain`   | –                      | UCSC liftOver chain file for coordinate remapping    |
| `--vep-forks`     | `4`                    | Number of parallel VEP forks                         |
| `--vep-custom`    | –                      | Passed to VEP `--custom`                             |
| `--vep-plugins`   | –                      | Passed to VEP `--plugin`                             |

---

## Architecture

``` text
vcf2maf/
├── constants.py   # VEP consequence priority table, MAF column list,
│                  # biotype rankings, VEP→MAF effect map
├── vcf2maf.py     # VCF → MAF  (VEP runner + annotation parser + MAF writer)
├── maf2vcf.py     # MAF → per-TN-pair VCF  (indel normalisation, samtools)
├── maf2maf.py     # MAF → MAF  (orchestrates maf2vcf + vcf2maf + merge)
├── vcf2vcf.py     # VCF → normalised VCF  (multiallelic split, liftOver)
└── pyproject.toml
```

---

## Differences from the Perl original

| Area             | Perl                           | Python                                |
| ---------------- | ------------------------------ | ------------------------------------- |
| Language         | Perl 5                         | Python ≥ 3.8                          |
| Dependencies     | CPAN modules                   | stdlib only                           |
| VEP interaction  | `system()` call                | `subprocess.run()`                    |
| Reference lookup | inline `samtools faidx`        | same via subprocess                   |
| Liftover         | `liftOver` binary              | same via subprocess                   |
| Parallelism      | `--vep-forks` forwarded to VEP | same                                  |
| ExAC columns     | present                        | preserved (empty; use gnomAD columns) |

Functional behaviour is identical for all standard use-cases. The Python code
explicitly documents each decision point that mirrors the Perl logic.

---

## Citation

If you use this tool in published research, please cite the original:

> Cyriac Kandoth. mskcc/vcf2maf: vcf2maf v1.6. (2020). doi:10.5281/zenodo.593251

## License

Apache-2.0 – same as the upstream mskcc/vcf2maf project.
