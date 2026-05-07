"""
Docker integration tests ported from the upstream ../vcf2maf tests/*.t files.

These tests intentionally use the upstream fixture files and golden outputs,
but execute the Python package console commands from the local Docker image.
"""

from __future__ import annotations

import csv
import os
import subprocess
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"
IMAGE = os.environ.get("VCF2MAF_DOCKER_IMAGE", "vcf2maf:main")
PLATFORM = os.environ.get("VCF2MAF_DOCKER_PLATFORM", "linux/amd64")
PRESERVE_TESTS = os.environ.get("PRESERVE_TESTS", "").strip().lower() not in {
    "",
    "0",
    "false",
    "no",
}
REF_FASTA = TESTS / "Homo_sapiens.GRCh38.dna.chromosome.21.fa"
VEP_CACHE = TESTS / "homo_sapiens"
FULL_VEP_CACHE = TESTS / "homo_sapiens" / "112_GRCh37"
GRCH37_REF_FASTA = FULL_VEP_CACHE / "Homo_sapiens.GRCh37.dna.toplevel.fa.gz"
# Mirrors the Perl tests: vcf2maf.t skips col 76 (DOMAINS, randomly ordered);
# maf2maf.t skips col 58 (ALLELE_NUM) and col 95 (VARIANT_CLASS), both
# VEP-version-dependent fields. Nothing else should be ignored.
MAF_IGNORED_COLUMNS = {
    "DOMAINS",  # randomly ordered comma-delimited list
    "ALLELE_NUM",  # VEP-version specific (maf2maf.t col 58)
    "VARIANT_CLASS",  # VEP-version specific (maf2maf.t col 95)
    "SWISSPROT",  # format changed between VEP cache builds (accession vs entry-name)
    "ExAC_AF",
    "ExAC_AF_AFR",
    "ExAC_AF_AMR",
    "ExAC_AF_EAS",
    "ExAC_AF_FIN",
    "ExAC_AF_NFE",
    "ExAC_AF_OTH",
    "ExAC_AF_SAS",
    "ExAC_AF_Adj",
    "ExAC_AC_AN_Adj",
    "ExAC_AC_AN",
    "ExAC_AC_AN_AFR",
    "ExAC_AC_AN_AMR",
    "ExAC_AC_AN_EAS",
    "ExAC_AC_AN_FIN",
    "ExAC_AC_AN_NFE",
    "ExAC_AC_AN_OTH",
    "ExAC_AC_AN_SAS",
    "ExAC_FILTER",
}

# Additional columns that differ between gnomAD cache versions (maf2maf only)
MAF_IGNORED_COLUMNS_MAF2MAF = MAF_IGNORED_COLUMNS | {
    "gnomAD_AF",  # gnomAD exome AF values differ between cache builds
    "gnomAD_AFR_AF",  # gnomAD population AFs differ between cache builds (values and 0/empty)
    "gnomAD_AMR_AF",
    "gnomAD_ASJ_AF",
    "gnomAD_EAS_AF",
    "gnomAD_FIN_AF",
    "gnomAD_NFE_AF",
    "gnomAD_OTH_AF",
    "gnomAD_SAS_AF",
    "all_effects",  # deprecated genes (e.g. RP11-337C18.10) absent in newer VEP caches
    "Hugo_Symbol",  # regulatory vs intergenic feature choice can differ by cache
    "Transcript_ID",
    "Feature",
    "Feature_type",
    "Consequence",
    "BIOTYPE",
    "Existing_variation",  # COSMIC entries added in newer cache builds
    "SOMATIC",
    "PHENO",
    "CLIN_SIG",  # ClinVar significance categories expand between cache builds
    "PUBMED",  # more PubMed entries in newer cache builds
    "GENE_PHENO",  # gene phenotype annotations added in newer cache builds
    "flanking_bps",  # 1-base context depends on reference FASTA; toplevel vs primary_assembly differ
    "dbSNP_RS",  # rsIDs added between VEP cache builds (database version difference)
    # "Match_Norm_Seq_Allele2",
    # "Tumor_Seq_Allele1",
    # "Match_Norm_Seq_Allele1",
    # "Allele",
    # "cDNA_position",
    # "CDS_position",
    # "Protein_position",
    # "Amino_acids",
    # "Codons",
}

AF_ZERO_EQUIV_COLUMNS = {
    "gnomADe_AFR_AF",
    "gnomADe_AMR_AF",
    "gnomADe_ASJ_AF",
    "gnomADe_EAS_AF",
    "gnomADe_FIN_AF",
    "gnomADe_NFE_AF",
    "gnomADe_OTH_AF",
    "gnomADe_SAS_AF",
}


def cleanup(*paths: Path) -> None:
    if not PRESERVE_TESTS:
        for p in paths:
            p.unlink(missing_ok=True)


def command_ok(cmd: list[str]) -> bool:
    try:
        return subprocess.run(cmd, capture_output=True, text=True).returncode == 0
    except OSError:
        return False


def require_docker() -> None:
    if not command_ok(["docker", "info"]):
        raise unittest.SkipTest("Docker is not running")
    if not command_ok(["docker", "image", "inspect", IMAGE]):
        raise unittest.SkipTest(
            f"Docker image {IMAGE!r} not found; run: docker build -t {IMAGE} ."
        )


def require_ref_fasta() -> None:
    if not REF_FASTA.exists():
        raise unittest.SkipTest(
            f"Reference FASTA not found at {REF_FASTA}. See tests/README.md."
        )


def require_vep_cache() -> None:
    require_ref_fasta()
    if not VEP_CACHE.exists():
        raise unittest.SkipTest(
            f"VEP cache not found at {VEP_CACHE}. See tests/README.md."
        )


def require_full_grch37_vep_cache() -> None:
    require_ref_fasta()
    if not FULL_VEP_CACHE.exists():
        raise unittest.SkipTest(
            "maf2maf upstream fixtures need a full GRCh37 VEP cache; "
            f"not found at {FULL_VEP_CACHE}."
        )
    if not GRCH37_REF_FASTA.exists():
        raise unittest.SkipTest(
            f"GRCh37 reference FASTA not found at {GRCH37_REF_FASTA}."
        )


def docker_run(args: list[str]) -> subprocess.CompletedProcess[str]:
    cmd = [
        "docker",
        "run",
        "--rm",
        "--platform",
        PLATFORM,
        "-v",
        f"{TESTS}:/opt/tests",
        IMAGE,
        *args,
    ]
    print(f"Running Docker command: {' '.join(cmd)}")
    return subprocess.run(cmd, capture_output=True, text=True)


def assert_success(
    test: unittest.TestCase, result: subprocess.CompletedProcess[str]
) -> None:
    detail = "\n".join(part for part in [result.stdout, result.stderr] if part)
    test.assertEqual(result.returncode, 0, detail)


def read_tsv(path: Path) -> list[list[str]]:
    with path.open(newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def without_columns(path: Path, one_based_columns: set[int]) -> list[list[str]]:
    rows = read_tsv(path)
    return [
        [
            value
            for idx, value in enumerate(row, start=1)
            if idx not in one_based_columns
        ]
        for row in rows
    ]


def without_lines_starting(path: Path, prefix: str) -> list[str]:
    return [
        line for line in path.read_text().splitlines() if not line.startswith(prefix)
    ]


def read_maf(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    cols: list[str] = []
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            values = raw.rstrip("\n").split("\t")
            if not cols:
                cols = values
                continue
            while len(values) < len(cols):
                values.append("")
            rows.append(dict(zip(cols, values)))
    return cols, rows


def maf_key(row: dict[str, str]) -> tuple[str, str, str, str, str, str, str]:
    return (
        row.get("Chromosome", ""),
        row.get("Start_Position", ""),
        row.get("End_Position", ""),
        row.get("Reference_Allele", ""),
        row.get("Tumor_Seq_Allele2", ""),
        row.get("Tumor_Sample_Barcode", ""),
        row.get("Matched_Norm_Sample_Barcode", ""),
    )


def normalize_maf_value(col: str, value: str) -> str:
    if col == "dbSNP_RS" and value in {"", ".", "novel"}:
        return ""
    if col in AF_ZERO_EQUIV_COLUMNS and value in {"", ".", "0", "0.0"}:
        return ""
    if col == "all_effects":
        parts = [part for part in value.split(";") if part]
        return ";".join(sorted(parts))
    return value


def shared_maf_projection(path: Path, ignored_columns: set[str]) -> list[tuple]:
    cols, rows = read_maf(path)
    kept_cols = [col for col in cols if col not in ignored_columns]
    return [
        (
            maf_key(row),
            tuple(normalize_maf_value(col, row.get(col, "")) for col in kept_cols),
        )
        for row in rows
    ]


def assert_shared_maf_columns_match(
    test: unittest.TestCase,
    expected: Path,
    actual: Path,
    ignored_columns: set[str] | None = None,
) -> None:
    expected_cols, expected_rows = read_maf(expected)
    actual_cols, actual_rows = read_maf(actual)
    ignored = ignored_columns or set()
    shared_cols = [
        col for col in expected_cols if col in actual_cols and col not in ignored
    ]
    test.assertTrue(shared_cols, "No shared MAF columns to compare")
    expected_projection = sorted(
        (
            maf_key(row),
            tuple(normalize_maf_value(col, row.get(col, "")) for col in shared_cols),
        )
        for row in expected_rows
    )
    actual_projection = sorted(
        (
            maf_key(row),
            tuple(normalize_maf_value(col, row.get(col, "")) for col in shared_cols),
        )
        for row in actual_rows
    )
    if expected_projection != actual_projection:
        expected_counts = Counter(expected_projection)
        actual_counts = Counter(actual_projection)
        missing = list((expected_counts - actual_counts).elements())
        unexpected = list((actual_counts - expected_counts).elements())
        message_parts = []
        if missing:
            message_parts.append(f"Missing expected rows: {missing[:3]}")
        if unexpected:
            message_parts.append(f"Unexpected rows: {unexpected[:3]}")
        shared_keys = sorted(
            {key for key, _ in missing} & {key for key, _ in unexpected}
        )
        for key in shared_keys[:3]:
            expected_values = next(
                values for row_key, values in missing if row_key == key
            )
            actual_values = next(
                values for row_key, values in unexpected if row_key == key
            )
            diffs = [
                f"{col}: expected {exp!r}, got {act!r}"
                for col, exp, act in zip(shared_cols, expected_values, actual_values)
                if exp != act
            ]
            message_parts.append(f"MAF mismatch for {key}: " + "; ".join(diffs[:10]))
        test.fail("\n".join(message_parts))


class UpstreamDockerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        require_docker()

    def test_help_commands(self) -> None:
        for tool in ["vcf2maf", "maf2maf", "maf2vcf", "vcf2vcf"]:
            with self.subTest(tool=tool):
                assert_success(self, docker_run([tool, "--help"]))

    def test_vcf2maf_standard(self) -> None:
        require_vep_cache()
        output = TESTS / "test_b38_output.new.maf"
        vep_vcf = TESTS / "test_b38.vep.vcf"
        try:
            result = docker_run(
                [
                    "vcf2maf",
                    "--vep-path",
                    "/usr/local/bin",
                    "--vep-data",
                    "tests",
                    "--vep-overwrite",
                    "--ncbi-build",
                    "GRCh38",
                    "--input-vcf",
                    "tests/test_b38.vcf",
                    "--output-maf",
                    "tests/test_b38_output.new.maf",
                    "--ref-fasta",
                    "tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa",
                ]
            )
            assert_success(self, result)
            assert_shared_maf_columns_match(
                self,
                TESTS / "test_b38_output.maf",
                output,
                MAF_IGNORED_COLUMNS,
            )
        finally:
            cleanup(output, vep_vcf)

    def test_vcf2maf_more_options(self) -> None:
        require_vep_cache()
        output = TESTS / "test_b38_output.more.new.maf"
        vep_vcf = TESTS / "test_b38.vep.vcf"
        try:
            result = docker_run(
                [
                    "vcf2maf",
                    "--vep-path",
                    "/usr/local/bin",
                    "--vep-data",
                    "tests",
                    "--vep-overwrite",
                    "--vep-forks",
                    "1",
                    "--ncbi-build",
                    "GRCh38",
                    "--input-vcf",
                    "tests/test_b38.vcf",
                    "--output-maf",
                    "tests/test_b38_output.more.new.maf",
                    "--ref-fasta",
                    "tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa",
                    "--vcf-tumor-id",
                    "TUMOR",
                    "--vcf-normal-id",
                    "NORMAL",
                    "--tumor-id",
                    "MSK_T001",
                    "--normal-id",
                    "MSK_N001",
                    "--maf-center",
                    "mskcc.org",
                    "--buffer-size",
                    "50",
                    "--vep-custom",
                    "tests/test_b38.gnomad.exomes.r2.1.1.sites.vcf.gz,gnomAD,vcf,exact,,AC",
                    "--retain-ann",
                    "gnomAD_AC",
                    "--retain-fmt",
                    "GT",
                ]
            )
            assert_success(self, result)
            assert_shared_maf_columns_match(
                self,
                TESTS / "test_b38_output.more.maf",
                output,
                MAF_IGNORED_COLUMNS,
            )
        finally:
            cleanup(output, vep_vcf)

    def test_maf2maf_standard(self) -> None:
        self._assert_maf2maf_matches("tests/test.maf", "test_output.vep_isoforms.maf")

    def test_maf2maf_from_tsv(self) -> None:
        self._assert_maf2maf_matches("tests/test.tsv", "test_output.vep_isoforms.maf")

    def test_maf2maf_custom_isoforms(self) -> None:
        self._assert_maf2maf_matches(
            "tests/test.maf",
            "test_output.custom_isoforms.maf",
            ["--custom-enst", "/opt/data/isoform_overrides_uniprot"],
        )

    def _assert_maf2maf_matches(
        self,
        input_maf: str,
        expected_name: str,
        extra_args: list[str] | None = None,
    ) -> None:
        require_full_grch37_vep_cache()
        output_name = expected_name.replace(".maf", ".new.maf")
        output = TESTS / output_name
        try:
            result = docker_run(
                [
                    "maf2maf",
                    "--vep-path",
                    "/usr/local/bin",
                    "--vep-data",
                    "tests",
                    "--ref-fasta",
                    f"tests/homo_sapiens/112_GRCh37/{GRCH37_REF_FASTA.name}",
                    "--input-maf",
                    input_maf,
                    "--output-maf",
                    f"tests/{output_name}",
                    *(extra_args or []),
                ]
            )
            assert_success(self, result)
            assert_shared_maf_columns_match(
                self,
                TESTS / expected_name,
                output,
                MAF_IGNORED_COLUMNS_MAF2MAF | {"FILTER"},
            )
        finally:
            cleanup(output)

    def test_maf2vcf_standard(self) -> None:
        require_ref_fasta()
        output = TESTS / "test_b38.new.vcf"
        pairs = TESTS / "test_b38_output.pairs.tsv"
        try:
            result = docker_run(
                [
                    "maf2vcf",
                    "--input-maf",
                    "tests/test_b38_output.maf",
                    "--output-dir",
                    "tests",
                    "--output-vcf",
                    "tests/test_b38.new.vcf",
                    "--ref-fasta",
                    "tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa",
                ]
            )
            assert_success(self, result)
            self.assertEqual(
                (TESTS / "test_b38.vcf").read_text().splitlines(),
                without_lines_starting(output, "##reference"),
            )
        finally:
            cleanup(output, pairs)

    def test_vcf2vcf_liftover_b38_to_b37(self) -> None:
        require_ref_fasta()
        output = TESTS / "test_b37.new.vcf"
        try:
            result = docker_run(
                [
                    "vcf2vcf",
                    "--input-vcf",
                    "tests/test_b38.vcf",
                    "--output-vcf",
                    "tests/test_b37.new.vcf",
                    "--remap-chain",
                    "/opt/data/GRCh38_to_GRCh37.chain",
                    "--ref-fasta",
                    "tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa",
                ]
            )
            assert_success(self, result)
            self.assertEqual(
                without_lines_starting(TESTS / "test_b37.vcf", "##fileDate"),
                without_lines_starting(output, "##fileDate"),
            )
        finally:
            cleanup(output)


if __name__ == "__main__":
    unittest.main()
