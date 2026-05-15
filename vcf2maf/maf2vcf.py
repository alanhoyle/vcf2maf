#!/usr/bin/env python3
"""
maf2vcf.py  –  Convert a MAF (or MAF-like TSV) into per-tumor/normal-pair VCFs,
suitable for re-annotation by vcf2maf.py.

Mirrors the logic of mskcc/vcf2maf  maf2vcf.pl.

Usage:
    python maf2vcf.py --input-maf INPUT.maf --output-dir vcfs/ --ref-fasta hg19.fa
"""

import argparse
import logging
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

log = logging.getLogger("maf2vcf")

_CHROM_ORDER = {"X": 23, "Y": 24, "MT": 25, "M": 25}


def _chrom_sort_key(chrom: str) -> tuple:
    c = chrom.upper().removeprefix("CHR")
    try:
        return (0, int(c))
    except ValueError:
        return (1, _CHROM_ORDER.get(c, 99))


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Convert a MAF into per-tumor/normal-pair VCFs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", default="",
                   help="Config file to load defaults from (default: ~/.vcf2maf.cfg)")
    p.add_argument("--input-maf", required=True, help="Input MAF file")
    p.add_argument(
        "--output-dir", required=True, help="Output directory for per-pair VCFs"
    )
    p.add_argument(
        "--output-vcf", default="", help="Output path for a single merged VCF"
    )
    p.add_argument(
        "--ref-fasta",
        default=os.path.expanduser(
            "~/.vep/homo_sapiens/112_GRCh38/"
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        ),
        help="Reference FASTA (must be samtools-indexed)",
    )
    p.add_argument(
        "--per-tn-vcfs",
        action="store_true",
        help="Write one VCF per tumor/normal pair (default: one merged VCF)",
    )
    p.add_argument(
        "--tumor-depth-col", default="t_depth", help="Column for tumor depth"
    )
    p.add_argument(
        "--tumor-vad-col", default="t_alt_count", help="Column for tumor alt depth"
    )
    p.add_argument(
        "--normal-depth-col", default="n_depth", help="Column for normal depth"
    )
    p.add_argument(
        "--normal-vad-col", default="n_alt_count", help="Column for normal alt depth"
    )
    p.add_argument("--samtools", default="samtools", help="Path to samtools binary")
    p.add_argument("--ncbi-build", default="GRCh38")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--debug", action="store_true", help="Set log level to DEBUG")
    return p


# ---------------------------------------------------------------------------
# MAF parsing helpers
# ---------------------------------------------------------------------------

REQUIRED_MAF_COLS = [
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Sample_Barcode",
]


def read_maf_header(path: str) -> Tuple[List[str], int]:
    """
    Return (column_names, header_line_number).
    Skips comment lines (starting with '#') and finds the first tab-separated
    line whose first field is one of the recognised MAF header fields.
    """
    header_fields = {
        "hugo_symbol",
        "chromosome",
        "tumor_sample_barcode",
        "start_position",
        "feature",
        "hugo_symbol",
    }
    with open(path) as fh:
        for lineno, raw in enumerate(fh, 1):
            if raw.startswith("#"):
                continue
            cols = raw.rstrip("\n\r").split("\t")
            if cols[0].lower() in header_fields or any(
                c.lower() in header_fields for c in cols
            ):
                return [c.strip() for c in cols], lineno
    raise ValueError(f"Could not find a valid MAF header in {path}")


def parse_maf_row(row: List[str], col_idx: Dict[str, int]) -> Dict[str, str]:
    """Map column names → values for one MAF row."""
    d: Dict[str, str] = {}
    for name, idx in col_idx.items():
        d[name] = row[idx] if idx < len(row) else ""
    return d


# ---------------------------------------------------------------------------
# Indel normalisation helpers
# ---------------------------------------------------------------------------


def maf_alleles_to_vcf(
    chrom: str, pos: int, ref: str, alt: str, flanking: str = ""
) -> Tuple[int, str, str]:
    """
    Convert MAF-style alleles (with '-' for gap) back to VCF-style left-aligned
    alleles.  Returns (vcf_pos, vcf_ref, vcf_alt).

    `flanking` should be the base immediately 5' of the variant in the reference
    (required only for indels where VCF needs an anchor base).
    """
    anchor = flanking[-1] if flanking else "N"  # preceding base as anchor

    if ref == "-":
        # Insertion: MAF pos is left of insertion; VCF pos is the anchor base
        return pos, anchor, anchor + alt

    if alt == "-":
        # Deletion: MAF start is first deleted base; VCF starts one base earlier
        return pos - 1, anchor + ref, anchor

    # Non-dash indel where lengths differ and no shared leading base: add anchor.
    # Mirrors the Perl maf2vcf.pl condition:
    #   $ref_len ne $al2_len and substr($ref,0,1) ne substr($al2,0,1)
    if len(ref) != len(alt) and ref[0] != alt[0]:
        return pos - 1, anchor + ref, anchor + alt

    # Pure substitution (same length, or shared leading base)
    return pos, ref, alt


def fetch_flanking_base(
    chrom: str, pos: int, ref_fasta: str, samtools: str = "samtools"
) -> str:
    """
    Use samtools faidx to retrieve the base immediately 5' of pos.
    Returns a single character, or 'N' on failure.
    """
    region = f"{chrom}:{pos-1}-{pos-1}"
    try:
        result = subprocess.run(
            [samtools, "faidx", ref_fasta, region],
            capture_output=True,
            text=True,
            check=True,
        )
        lines = result.stdout.strip().split("\n")
        if len(lines) >= 2:
            return lines[1].strip().upper() or "N"
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return "N"


def fetch_reference_base(
    chrom: str, pos: int, ref_fasta: str, samtools: str = "samtools"
) -> str:
    region = f"{chrom}:{pos}-{pos}"
    try:
        result = subprocess.run(
            [samtools, "faidx", ref_fasta, region],
            capture_output=True,
            text=True,
            check=True,
        )
        lines = result.stdout.strip().split("\n")
        if len(lines) >= 2:
            return lines[1].strip().upper() or "N"
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return "N"


# ---------------------------------------------------------------------------
# VCF writer
# ---------------------------------------------------------------------------

VCF_HEADER_TEMPLATE = """\
##fileformat=VCFv4.2
##reference={ref_fasta}
##contig=<ID={chrom},length={chrom_length}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths of REF and ALT(s) in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL
"""


def chrom_length_from_fai(ref_fasta: str, chrom: str) -> str:
    fai = ref_fasta + ".fai"
    try:
        with open(fai) as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2 and fields[0] == chrom:
                    return fields[1]
    except OSError:
        pass
    return "."


def write_vcf(
    path: str,
    variants: List[Dict],
    tumor_id: str,
    normal_id: str,
    ref_fasta: str,
    assembly: str,
) -> None:
    """Write a VCF file from a list of variant dicts."""
    with open(path, "w") as fh:
        header = VCF_HEADER_TEMPLATE.format(
            ref_fasta=ref_fasta,
            chrom=variants[0]["vcf_chrom"] if variants else "*",
            chrom_length=chrom_length_from_fai(
                ref_fasta, variants[0]["vcf_chrom"] if variants else "*"
            ),
        )
        # Replace TUMOR/NORMAL placeholder with real IDs
        header = header.replace("\tTUMOR\tNORMAL", f"\t{tumor_id}\t{normal_id}")
        fh.write(header)

        for v in variants:
            chrom = v["vcf_chrom"]
            pos = v["vcf_pos"]
            ref = v["vcf_ref"]
            alt = v["vcf_alt"]
            filt = v.get("filter", ".") or "."
            vcf_id = v.get("variant_id", ".") or "."

            info_str = "."

            fmt_fields = "GT:AD:DP"
            tum_sample = v.get("tumor_format")
            nrm_sample = v.get("normal_format")
            if not tum_sample or not nrm_sample:

                def fmt_gt(
                    depth: str, ref_count: str, alt_count: str, is_alt: bool
                ) -> str:
                    gt = "0/1" if is_alt else "0/0"
                    ref_count = ref_count or "."
                    alt_count = alt_count or "."
                    depth = depth or "."
                    return f"{gt}:{ref_count},{alt_count}:{depth}"

                tum_sample = fmt_gt(
                    v.get("t_depth", ""),
                    v.get("t_ref_count", ""),
                    v.get("t_alt_count", ""),
                    True,
                )
                nrm_sample = fmt_gt(
                    v.get("n_depth", ""),
                    v.get("n_ref_count", ""),
                    v.get("n_alt_count", ""),
                    False,
                )

            fh.write(
                f"{chrom}\t{pos}\t{vcf_id}\t{ref}\t{alt}\t.\t{filt}\t"
                f"{info_str}\t{fmt_fields}\t{tum_sample}\t{nrm_sample}\n"
            )


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------


def maf2vcf(args: argparse.Namespace) -> None:
    """Convert a MAF to per-TN-pair VCFs."""

    if not os.path.isfile(args.input_maf):
        sys.exit(f"ERROR: --input-maf not found: {args.input_maf}")
    if not os.path.isfile(args.ref_fasta):
        sys.exit(f"ERROR: --ref-fasta not found: {args.ref_fasta}")

    log.info("maf2vcf: %s → %s", args.input_maf, args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    col_names, header_lineno = read_maf_header(args.input_maf)
    col_idx = {c.lower(): i for i, c in enumerate(col_names)}
    log.info("MAF header found at line %d: %d columns", header_lineno, len(col_names))

    # Validate required columns
    for req in REQUIRED_MAF_COLS:
        if req.lower() not in col_idx:
            sys.exit(f"ERROR: Required MAF column '{req}' not found in header")

    # Collect all variants, grouped by (tumor_id, normal_id)
    tn_variants: Dict[Tuple[str, str], List[Dict]] = defaultdict(list)
    mismatch_rows: List[str] = []

    with open(args.input_maf) as fh:
        for lineno, raw in enumerate(fh, 1):
            if lineno <= header_lineno:
                continue
            if raw.startswith("#"):
                continue
            line = raw.rstrip("\n\r")
            if not line:
                continue

            cols_raw = line.split("\t")
            row = parse_maf_row(cols_raw, col_idx)

            chrom = row.get("chromosome", "")
            pos_str = row.get("start_position", "")
            ref = re.sub(r"^[\?0]+$", "", row.get("reference_allele", ""))
            al1 = row.get("tumor_seq_allele1", "")
            al2 = row.get("tumor_seq_allele2", "")
            n_al1 = row.get("match_norm_seq_allele1", "")
            n_al2 = row.get("match_norm_seq_allele2", "")
            tumor_id = row.get("tumor_sample_barcode", "TUMOR")
            normal_id = row.get(
                "matched_norm_sample_barcode",
                row.get("normal_sample_barcode", "NORMAL"),
            )

            if not pos_str.isdigit():
                continue
            pos = int(pos_str)

            al1 = al1 or ref
            al2 = al2 or ref
            n_al1 = n_al1 or ref
            n_al2 = n_al2 or ref

            # Upstream maf2vcf.pl blanks dash-like indel alleles before VCF
            # anchoring, then restores proper VCF alleles with the prefix base.
            ref, al1, al2, n_al1, n_al2 = [
                re.sub(r"^[\?\-0]+$", "", allele)
                for allele in (ref, al1, al2, n_al1, n_al2)
            ]

            if al2 == ref:
                al1, al2 = al2, al1
            if n_al2 == ref:
                n_al1, n_al2 = n_al2, n_al1

            alt = al2
            if alt == ref:
                log.debug("Skip %s:%d — alt equals ref after normalization", chrom, pos)
                continue

            ref_len, al1_len, al2_len = len(ref), len(al1), len(al2)
            needs_anchor = (
                ref_len == 0
                or al1_len == 0
                or al2_len == 0
                or (ref_len != al2_len and ref[0] != al2[0])
            )
            vcf_pos = pos
            if needs_anchor:
                if ref == "":
                    prefix = fetch_reference_base(
                        chrom, pos, args.ref_fasta, args.samtools
                    )
                else:
                    prefix = fetch_flanking_base(
                        chrom, pos, args.ref_fasta, args.samtools
                    )
                    vcf_pos -= 1
                log.debug("Fetched anchor base %r for %s:%d", prefix, chrom, pos)
                ref, al1, al2, n_al1, n_al2 = [
                    prefix + allele for allele in (ref, al1, al2, n_al1, n_al2)
                ]

            alleles: List[str] = []
            allele_indices: Dict[str, int] = {}
            for allele in (ref, al2, al1, n_al2, n_al1):
                if allele not in allele_indices:
                    allele_indices[allele] = len(alleles)
                    alleles.append(allele)

            vcf_ref = alleles[0]
            vcf_alt = ",".join(alleles[1:])

            t_gt = "0/1"
            if allele_indices.get(al1) != 0:
                t_gt = f"{allele_indices[al2]}/{allele_indices[al1]}"

            n_gt = "0/0"
            if allele_indices.get(n_al1) != 0:
                n_gt = f"{allele_indices[n_al2]}/{allele_indices[n_al1]}"
            if allele_indices.get(n_al2) != 0:
                n_gt = f"{allele_indices[n_al1]}/{allele_indices[n_al2]}"

            t_depth = row.get(args.tumor_depth_col.lower(), "")
            t_ref_count = row.get("t_ref_count", "")
            t_alt_count = row.get(args.tumor_vad_col.lower(), "")
            n_depth = row.get(args.normal_depth_col.lower(), "")
            n_ref_count = row.get("n_ref_count", "")
            n_alt_count = row.get(args.normal_vad_col.lower(), "")

            extra_alt_depths = ["." for _ in alleles[2:]]
            t_ad = ",".join([t_ref_count or ".", t_alt_count or ".", *extra_alt_depths])
            n_ad = ",".join([n_ref_count or ".", n_alt_count or ".", *extra_alt_depths])

            v: Dict = dict(row)
            v.update(
                {
                    "vcf_chrom": chrom,
                    "vcf_pos": vcf_pos,
                    "vcf_ref": vcf_ref,
                    "vcf_alt": vcf_alt,
                    "t_depth": t_depth,
                    "t_alt_count": t_alt_count,
                    "t_ref_count": t_ref_count,
                    "n_depth": n_depth,
                    "n_alt_count": n_alt_count,
                    "n_ref_count": n_ref_count,
                    "tumor_format": f"{t_gt}:{t_ad}:{t_depth or '.'}",
                    "normal_format": f"{n_gt}:{n_ad}:{n_depth or '.'}",
                }
            )

            tn_variants[(tumor_id, normal_id)].append(v)

    total_variants = sum(len(v) for v in tn_variants.values())
    log.info(
        "Collected %d variant(s) across %d TN pair(s)",
        total_variants, len(tn_variants),
    )

    # Write VCFs
    written_vcfs: List[str] = []
    for (tumor_id, normal_id), variants in tn_variants.items():
        if args.per_tn_vcfs:
            # The Perl implementation writes every per-TN input row. Some
            # maf2maf fixtures intentionally have repeated loci with different
            # normal alleles/counts, and deduping loses those distinctions.
            unique_variants = list(variants)
        else:
            seen_vars: set = set()
            unique_variants = []
            for v in variants:
                key = (v["vcf_chrom"], v["vcf_pos"], v["vcf_ref"], v["vcf_alt"])
                if key not in seen_vars:
                    seen_vars.add(key)
                    unique_variants.append(v)

        # Sort by chrom (numerically for 1-22, then X/Y/MT) / pos
        unique_variants.sort(
            key=lambda v: (_chrom_sort_key(v["vcf_chrom"]), v["vcf_pos"])
        )

        safe_tumor = re.sub(r"[^\w\-]", "_", tumor_id)
        safe_normal = re.sub(r"[^\w\-]", "_", normal_id)
        if not args.per_tn_vcfs and args.output_vcf:
            vcf_path = args.output_vcf
        else:
            vcf_path = os.path.join(
                args.output_dir, f"{safe_tumor}_vs_{safe_normal}.vcf"
            )
        write_vcf(
            vcf_path,
            unique_variants,
            tumor_id,
            normal_id,
            args.ref_fasta,
            args.ncbi_build,
        )
        written_vcfs.append(vcf_path)
        log.info("Wrote %d variants to %s", len(unique_variants), vcf_path)
        if not args.per_tn_vcfs:
            break

    return written_vcfs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default="")
    pre.add_argument("--debug", action="store_true")
    pre_args, _ = pre.parse_known_args()

    logging.basicConfig(
        level=logging.DEBUG if pre_args.debug else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    try:
        from .config import load_config
    except ImportError:
        from config import load_config  # type: ignore[no-redef]

    parser = build_parser()
    cfg = load_config(pre_args.config or None)
    if cfg:
        parser.set_defaults(**cfg)
    args = parser.parse_args()
    maf2vcf(args)


if __name__ == "__main__":
    main()
