#!/usr/bin/env python3
"""
maf2maf.py  –  Re-annotate a MAF by running maf2vcf then vcf2maf on each
tumor/normal pair, and merging the results back into a single MAF.

Mirrors mskcc/vcf2maf  maf2maf.pl.

Usage:
    python maf2maf.py --input-maf INPUT.maf --output-maf OUTPUT.vep.maf
"""

import argparse
import logging
import os
import re
import shutil
import sys
import tempfile
from importlib import import_module
from pathlib import Path
from typing import List

try:
    m2v = import_module(".maf2vcf", __package__)
    v2m = import_module(".vcf2maf", __package__)
except (ImportError, TypeError):
    m2v = import_module("maf2vcf")
    v2m = import_module("vcf2maf")

log = logging.getLogger("maf2maf")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Re-annotate a MAF via maf2vcf + vcf2maf.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", default="",
                   help="Config file to load defaults from (default: ~/.vcf2maf.cfg)")
    p.add_argument("--input-maf", required=True)
    p.add_argument("--output-maf", required=True)

    # Pass-through to maf2vcf
    p.add_argument(
        "--ref-fasta",
        default=os.path.expanduser(
            "~/.vep/homo_sapiens/112_GRCh38/"
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        ),
    )
    p.add_argument("--tumor-depth-col", default="t_depth")
    p.add_argument("--tumor-vad-col", default="t_alt_count")
    p.add_argument("--normal-depth-col", default="n_depth")
    p.add_argument("--normal-vad-col", default="n_alt_count")
    p.add_argument("--samtools", default="samtools")

    # Pass-through to vcf2maf
    p.add_argument("--vep-path", default=os.path.expanduser("~/miniconda3/bin"))
    p.add_argument("--vep-data", default=os.path.expanduser("~/.vep"))
    p.add_argument("--vep-forks", type=int, default=4)
    p.add_argument("--buffer-size", type=int, default=5000)
    p.add_argument("--cache-version", default="")
    p.add_argument("--vep-custom", default="")
    p.add_argument("--vep-config", default="")
    p.add_argument("--vep-plugins", default="")
    p.add_argument("--vep-overwrite", action="store_true")
    p.add_argument(
        "--vep-log-cmd",
        action="store_true",
        help="Log the VEP command in shell-style multi-line format (easier to copy/re-run).",
    )
    p.add_argument(
        "--vep-stats",
        nargs="?",
        const=True,
        default=None,
        help="Path for VEP summary stats file (HTML). "
        "Omit to suppress stats (--no_stats). "
        "Pass without a value to use VEP's default stats filename.",
    )
    p.add_argument(
        "--vep-stats-text",
        action="store_true",
        default=False,
        help="Also write a plain-text stats file (passes --stats_text to VEP).",
    )
    p.add_argument(
        "--vep-stats-html",
        action="store_true",
        default=False,
        help="Also write an HTML stats file (passes --stats_html to VEP).",
    )
    p.add_argument("--inhibit-vep", action="store_true")
    p.add_argument("--ncbi-build", default="GRCh38")
    p.add_argument("--species", default="homo_sapiens")
    p.add_argument("--maf-center", default=".")
    p.add_argument("--min-hom-vaf", type=float, default=0.7)
    p.add_argument("--max-subpop-af", type=float, default=0.0004)
    p.add_argument("--custom-enst", default="")
    p.add_argument("--retain-info", default="")
    p.add_argument("--retain-fmt", default="")
    p.add_argument("--retain-ann", default="")
    p.add_argument("--retain-ann-all", action="store_true", default=False,
                   help="Retain all VEP CSQ annotation fields as extra MAF columns")
    p.add_argument(
        "--retain-cols",
        default="",
        help="Comma-separated list of columns from the input MAF to carry over "
        "unchanged into the output MAF",
    )
    p.add_argument("--remap-chain", default="")
    p.add_argument("--liftover-exec", default="liftOver")
    p.add_argument("--verbose", action="store_true")
    return p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def maf2maf(args: argparse.Namespace) -> None:
    if not os.path.isfile(args.input_maf):
        sys.exit(f"ERROR: --input-maf not found: {args.input_maf}")

    log.info("maf2maf: %s → %s", args.input_maf, args.output_maf)
    tmp_dir = tempfile.mkdtemp(prefix="maf2maf_")
    log.info("Temporary directory: %s", tmp_dir)

    try:
        # ------------------------------------------------------------------
        # Step 1: maf2vcf – produce per-TN-pair VCFs
        # ------------------------------------------------------------------
        vcf_dir = os.path.join(tmp_dir, "vcfs")
        m2v_args = argparse.Namespace(
            input_maf=args.input_maf,
            output_dir=vcf_dir,
            ref_fasta=args.ref_fasta,
            per_tn_vcfs=True,
            tumor_depth_col=args.tumor_depth_col,
            tumor_vad_col=args.tumor_vad_col,
            normal_depth_col=args.normal_depth_col,
            normal_vad_col=args.normal_vad_col,
            samtools=args.samtools,
            ncbi_build=args.ncbi_build,
            verbose=args.verbose,
        )
        vcf_paths = m2v.maf2vcf(m2v_args)

        if not vcf_paths:
            sys.exit("ERROR: maf2vcf produced no VCF files")
        log.info("Step 1 complete: %d per-TN VCF(s) produced", len(vcf_paths))

        # ------------------------------------------------------------------
        # Step 2: vcf2maf – annotate each VCF
        # ------------------------------------------------------------------
        maf_parts: List[str] = []
        for vcf_path in vcf_paths:
            # Derive tumor/normal IDs from VCF filename: TUMOR_vs_NORMAL.vcf
            base = os.path.basename(vcf_path).replace(".vcf", "")
            parts = base.split("_vs_")
            tumor_id = parts[0] if parts else "TUMOR"
            normal_id = parts[1] if len(parts) > 1 else "NORMAL"

            part_maf = vcf_path.replace(".vcf", ".maf")

            v2m_args = argparse.Namespace(
                input_vcf=vcf_path,
                output_maf=part_maf,
                tumor_id=tumor_id,
                normal_id=normal_id,
                vcf_tumor_id=tumor_id,
                vcf_normal_id=normal_id,
                vep_path=args.vep_path,
                vep_data=args.vep_data,
                vep_forks=args.vep_forks,
                buffer_size=args.buffer_size,
                cache_version=args.cache_version,
                vep_custom=args.vep_custom,
                vep_config=args.vep_config,
                vep_plugins=args.vep_plugins,
                vep_overwrite=args.vep_overwrite,
                vep_log_cmd=args.vep_log_cmd,
                vep_stats=args.vep_stats,
                vep_stats_text=args.vep_stats_text,
                vep_stats_html=args.vep_stats_html,
                ref_fasta=args.ref_fasta,
                species=args.species,
                ncbi_build=args.ncbi_build,
                maf_center=args.maf_center,
                min_hom_vaf=args.min_hom_vaf,
                max_subpop_af=args.max_subpop_af,
                retain_info=args.retain_info,
                retain_fmt=args.retain_fmt,
                retain_ann=args.retain_ann,
                retain_ann_all=args.retain_ann_all,
                custom_enst=args.custom_enst,
                remap_chain=args.remap_chain,
                liftover_exec=args.liftover_exec,
                inhibit_vep=args.inhibit_vep,
                any_allele=False,
                online=False,
                verbose=args.verbose,
            )
            log.info(
                "Step 2 [%d/%d]: annotating %s (tumor=%s, normal=%s)",
                len(maf_parts) + 1, len(vcf_paths), vcf_path, tumor_id, normal_id,
            )
            v2m.vcf2maf(v2m_args)
            maf_parts.append(part_maf)

        # ------------------------------------------------------------------
        # Step 3: Merge per-TN MAFs → output MAF
        # ------------------------------------------------------------------
        log.info("Step 3: merging %d partial MAF(s) → %s", len(maf_parts), args.output_maf)
        _merge_mafs(maf_parts, args.output_maf, args.retain_cols, args.input_maf)
        log.info("Final MAF written to %s", args.output_maf)

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def _merge_mafs(
    maf_parts: List[str], output_maf: str, retain_cols: str, original_maf: str
) -> None:
    """
    Concatenate per-TN MAFs, writing the header only once.
    Optionally carry over user-specified columns from the original MAF
    (keyed on Chromosome + Start_Position + Reference_Allele + Tumor_Sample_Barcode).
    """
    retain = [c.strip() for c in retain_cols.split(",") if c.strip()]

    # Build lookup from original MAF for retain columns
    orig_lookup: dict = {}
    if retain:
        try:
            orig_col_names, orig_hdr_lineno = m2v.read_maf_header(original_maf)
            orig_col_idx = {c.lower(): i for i, c in enumerate(orig_col_names)}
            with open(original_maf) as fh:
                for lineno, raw in enumerate(fh, 1):
                    if lineno <= orig_hdr_lineno or raw.startswith("#"):
                        continue
                    parts = raw.rstrip("\n\r").split("\t")
                    chrom = _get(parts, orig_col_idx, "chromosome")
                    pos = _get(parts, orig_col_idx, "start_position")
                    ref = _get(parts, orig_col_idx, "reference_allele")
                    tsb = _get(parts, orig_col_idx, "tumor_sample_barcode")
                    key = (chrom, pos, ref, tsb)
                    orig_lookup[key] = {
                        c: _get(parts, orig_col_idx, c.lower()) for c in retain
                    }
        except Exception as exc:
            log.warning("Could not load retain columns from original MAF: %s", exc)

    header_written = False
    out_col_names: List[str] = []

    with open(output_maf, "w") as out_fh:
        for part in maf_parts:
            with open(part) as in_fh:
                for raw in in_fh:
                    if raw.startswith("#version"):
                        if not header_written:
                            out_fh.write(raw)
                        continue

                    if not raw.startswith("#") and not header_written:
                        # Header row
                        out_col_names = raw.rstrip("\n\r").split("\t")
                        # Append extra retain cols not already present
                        existing_set = set(c.lower() for c in out_col_names)
                        for rc in retain:
                            if rc.lower() not in existing_set:
                                out_col_names.append(rc)
                        out_fh.write("\t".join(out_col_names) + "\n")
                        header_written = True
                        continue

                    if raw.startswith("#"):
                        continue  # skip secondary header / comment lines

                    # Data row
                    row_parts = raw.rstrip("\n\r").split("\t")
                    # Pad short rows
                    while len(row_parts) < len(out_col_names):
                        row_parts.append("")

                    if retain and orig_lookup:
                        col_map = {c.lower(): i for i, c in enumerate(out_col_names)}
                        chrom = _get_by_map(row_parts, col_map, "chromosome")
                        pos = _get_by_map(row_parts, col_map, "start_position")
                        ref = _get_by_map(row_parts, col_map, "reference_allele")
                        tsb = _get_by_map(row_parts, col_map, "tumor_sample_barcode")
                        key = (chrom, pos, ref, tsb)
                        extra = orig_lookup.get(key, {})
                        for rc in retain:
                            rc_idx = col_map.get(rc.lower())
                            if rc_idx is not None:
                                row_parts[rc_idx] = extra.get(rc, row_parts[rc_idx])

                    out_fh.write("\t".join(row_parts) + "\n")


def _get(parts: list, idx: dict, key: str) -> str:
    i = idx.get(key)
    if i is not None and i < len(parts):
        return parts[i]
    return ""


def _get_by_map(parts: list, idx: dict, key: str) -> str:
    i = idx.get(key)
    if i is not None and i < len(parts):
        return parts[i]
    return ""


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    try:
        from .config import load_config
    except ImportError:
        from config import load_config  # type: ignore[no-redef]

    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", default="")
    pre_args, _ = pre.parse_known_args()

    parser = build_parser()
    cfg = load_config(pre_args.config or None)
    if cfg:
        parser.set_defaults(**cfg)
    args = parser.parse_args()
    maf2maf(args)


if __name__ == "__main__":
    main()
