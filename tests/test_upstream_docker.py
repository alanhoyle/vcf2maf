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
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"
IMAGE = os.environ.get("VCF2MAF_DOCKER_IMAGE", "vcf2maf:main")
PLATFORM = os.environ.get("VCF2MAF_DOCKER_PLATFORM", "linux/amd64")
REF_FASTA = TESTS / "Homo_sapiens.GRCh38.dna.chromosome.21.fa"


def command_ok(cmd: list[str]) -> bool:
    try:
        return subprocess.run(cmd, capture_output=True, text=True).returncode == 0
    except OSError:
        return False


def require_docker() -> None:
    if not command_ok(["docker", "info"]):
        raise unittest.SkipTest("Docker is not running")
    if not command_ok(["docker", "image", "inspect", IMAGE]):
        raise unittest.SkipTest(f"Docker image {IMAGE!r} not found; run: docker build -t {IMAGE} .")


def require_ref_fasta() -> None:
    if not REF_FASTA.exists():
        raise unittest.SkipTest(
            f"Reference FASTA not found at {REF_FASTA}. See tests/README.md."
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
    return subprocess.run(cmd, capture_output=True, text=True)


def assert_success(test: unittest.TestCase, result: subprocess.CompletedProcess[str]) -> None:
    detail = "\n".join(part for part in [result.stdout, result.stderr] if part)
    test.assertEqual(result.returncode, 0, detail)


def read_tsv(path: Path) -> list[list[str]]:
    with path.open(newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def without_columns(path: Path, one_based_columns: set[int]) -> list[list[str]]:
    rows = read_tsv(path)
    return [
        [value for idx, value in enumerate(row, start=1) if idx not in one_based_columns]
        for row in rows
    ]


def without_lines_starting(path: Path, prefix: str) -> list[str]:
    return [line for line in path.read_text().splitlines() if not line.startswith(prefix)]


class UpstreamDockerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        require_docker()

    def test_help_commands(self) -> None:
        for tool in ["vcf2maf", "maf2maf", "maf2vcf", "vcf2vcf"]:
            with self.subTest(tool=tool):
                assert_success(self, docker_run([tool, "--help"]))

    def test_vcf2maf_standard(self) -> None:
        require_ref_fasta()
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
            self.assertEqual(
                without_columns(TESTS / "test_b38_output.maf", {76}),
                without_columns(output, {76}),
            )
        finally:
            output.unlink(missing_ok=True)
            vep_vcf.unlink(missing_ok=True)

    def test_vcf2maf_more_options(self) -> None:
        require_ref_fasta()
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
            self.assertEqual(
                without_columns(TESTS / "test_b38_output.more.maf", {76}),
                without_columns(output, {76}),
            )
        finally:
            output.unlink(missing_ok=True)
            vep_vcf.unlink(missing_ok=True)

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
        require_ref_fasta()
        output = TESTS / "test_output.vep_isoforms.new.maf"
        try:
            result = docker_run(
                [
                    "maf2maf",
                    "--vep-path",
                    "/usr/local/bin",
                    "--vep-data",
                    "tests",
                    "--ref-fasta",
                    "tests/Homo_sapiens.GRCh38.dna.chromosome.21.fa",
                    "--input-maf",
                    input_maf,
                    "--output-maf",
                    "tests/test_output.vep_isoforms.new.maf",
                    *(extra_args or []),
                ]
            )
            assert_success(self, result)
            self.assertEqual(
                without_columns(TESTS / expected_name, {58, 95}),
                without_columns(output, {58, 95}),
            )
        finally:
            output.unlink(missing_ok=True)

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
            output.unlink(missing_ok=True)
            pairs.unlink(missing_ok=True)

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
            output.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
