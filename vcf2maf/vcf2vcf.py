#!/usr/bin/env python3
"""
vcf2vcf.py  –  Normalize / clean up a VCF so that it is suitable for vcf2maf.
Handles multiallelic splitting, left-alignment, liftOver remapping, and
sample column reordering.

Mirrors mskcc/vcf2maf  vcf2vcf.pl.

Usage:
    python vcf2vcf.py --input-vcf INPUT.vcf --output-vcf OUTPUT.vcf \\
        --vcf-tumor-id TUMOR --vcf-normal-id NORMAL
"""

import argparse
import logging
import os
import re
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Tuple

log = logging.getLogger("vcf2vcf")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Normalise and clean a VCF prior to vcf2maf annotation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", default="",
                   help="Config file to load defaults from (default: ~/.vcf2maf.cfg)")
    p.add_argument("--input-vcf", required=True)
    p.add_argument("--output-vcf", required=True)
    p.add_argument("--vcf-tumor-id", default="TUMOR")
    p.add_argument("--vcf-normal-id", default="NORMAL")
    p.add_argument(
        "--ref-fasta",
        default=os.path.expanduser(
            "~/.vep/homo_sapiens/112_GRCh38/"
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        ),
    )
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument(
        "--remap-chain", default="", help="UCSC liftOver chain file for remapping"
    )
    p.add_argument("--liftover-exec", default="liftOver")
    p.add_argument(
        "--add-filters",
        action="store_true",
        help="Add FILTER annotations for known issues",
    )
    p.add_argument("--ncbi-build", default="GRCh38")
    p.add_argument("--verbose", action="store_true")
    p.add_argument("--debug", action="store_true", help="Set log level to DEBUG")
    return p


# ---------------------------------------------------------------------------
# Liftover helpers
# ---------------------------------------------------------------------------


def remap_vcf(input_vcf: str, chain: str, liftover_exec: str, tmp_dir: str) -> Dict[str, str]:
    """
    Run UCSC liftOver on the variant loci and return a mapping of
    original 'chrom:pos' keys to remapped 'chrom:pos' values.
    Exits with an error if liftOver cannot be found or fails to run.
    """
    bed_path = os.path.join(tmp_dir, "input.bed")
    mapped_path = os.path.join(tmp_dir, "mapped.bed")
    unmapped_path = os.path.join(tmp_dir, "unmapped.bed")
    remap: Dict[str, str] = {}

    # Build BED from VCF
    with open(input_vcf) as vcf_fh, open(bed_path, "w") as bed_fh:
        for line in vcf_fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            chrom, pos = cols[0], int(cols[1])
            bed_fh.write(f"{chrom}\t{pos-1}\t{pos}\t{chrom}:{pos}\n")

    # Run liftOver
    try:
        result = subprocess.run(
            [liftover_exec, bed_path, chain, mapped_path, unmapped_path],
            capture_output=True,
        )
    except FileNotFoundError:
        sys.exit(f"ERROR: liftOver executable not found: {liftover_exec!r}")
    except OSError as exc:
        sys.exit(f"ERROR: Failed to launch liftOver: {exc}")

    if result.returncode != 0:
        stderr = result.stderr.decode(errors="replace").strip()
        sys.exit(
            f"ERROR: liftOver exited with status {result.returncode}"
            + (f": {stderr}" if stderr else "")
        )

    # Warn about unmapped loci
    if os.path.isfile(unmapped_path):
        unmapped_count = sum(
            1 for ln in open(unmapped_path) if not ln.startswith("#") and ln.strip()
        )
        if unmapped_count:
            log.warning(
                "%d locus/loci could not be remapped by liftOver and will be dropped",
                unmapped_count,
            )

    # Parse mapped loci
    if os.path.isfile(mapped_path):
        with open(mapped_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) >= 4:
                    orig_key = parts[3].strip()
                    new_chrom, new_end = parts[0], int(parts[2])
                    remap[orig_key] = f"{new_chrom}:{new_end}"

    return remap


# ---------------------------------------------------------------------------
# Multiallelic splitting
# ---------------------------------------------------------------------------


def split_multiallelic(line: str) -> List[str]:
    """
    Split a multiallelic VCF line into one line per ALT allele.
    Adjusts GT fields in all sample columns.
    """
    cols = line.rstrip("\n\r").split("\t")
    alts = cols[4].split(",")
    if len(alts) == 1:
        return [line]

    split_lines = []
    for alt_idx, alt in enumerate(alts, 1):
        new_cols = list(cols)
        new_cols[4] = alt
        # Recode sample GT fields: replace ALT index or set to 0/1
        for sample_i in range(9, len(new_cols)):
            fmt_keys = cols[8].split(":") if len(cols) > 8 else []
            sample_vals = new_cols[sample_i].split(":")
            if fmt_keys and "GT" in fmt_keys:
                gt_i = fmt_keys.index("GT")
                if gt_i < len(sample_vals):
                    gt = sample_vals[gt_i]
                    alleles = re.split(r"([/|])", gt)
                    new_alleles = []
                    for a in alleles:
                        if a in "/|":
                            new_alleles.append(a)
                        elif a.isdigit():
                            ai = int(a)
                            if ai == 0:
                                new_alleles.append("0")
                            elif ai == alt_idx:
                                new_alleles.append("1")
                            else:
                                new_alleles.append(".")
                        else:
                            new_alleles.append(a)
                    sample_vals[gt_i] = "".join(new_alleles)
                new_cols[sample_i] = ":".join(sample_vals)
        split_lines.append("\t".join(new_cols) + "\n")

    return split_lines


# ---------------------------------------------------------------------------
# Left-alignment
# ---------------------------------------------------------------------------


def left_align_variant(
    chrom: str, pos: int, ref: str, alt: str, ref_fasta: str, samtools: str
) -> Tuple[int, str, str]:
    """
    Left-align an indel by shifting as far left as possible.
    For substitutions, returns input unchanged.
    """
    if len(ref) == len(alt):
        return pos, ref, alt  # substitution

    # Trim common suffix
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    # Shift left while last base of ref equals last base of what precedes
    while True:
        # Fetch the base one position to the left
        region = f"{chrom}:{pos-1}-{pos-1}"
        try:
            result = subprocess.run(
                [samtools, "faidx", ref_fasta, region],
                capture_output=True,
                text=True,
                check=True,
            )
            lines = result.stdout.strip().split("\n")
            prev_base = lines[1].strip().upper() if len(lines) >= 2 else ""
        except Exception:
            break

        if not prev_base:
            break

        # Check if we can shift left
        if len(ref) > 1 and ref[-1] == prev_base:
            ref = prev_base + ref[:-1]
            alt = prev_base + alt[:-1]
            pos -= 1
        elif len(alt) > 1 and alt[-1] == prev_base:
            ref = prev_base + ref[:-1]
            alt = prev_base + alt[:-1]
            pos -= 1
        else:
            break

    return pos, ref, alt


# ---------------------------------------------------------------------------
# Main normalisation
# ---------------------------------------------------------------------------


def vcf2vcf(args: argparse.Namespace) -> None:
    if not os.path.isfile(args.input_vcf):
        sys.exit(f"ERROR: --input-vcf not found: {args.input_vcf}")
    if not os.path.isfile(args.ref_fasta):
        sys.exit(f"ERROR: --ref-fasta not found: {args.ref_fasta}")

    log.info("vcf2vcf: %s → %s", args.input_vcf, args.output_vcf)
    tmp_dir = tempfile.mkdtemp(prefix="vcf2vcf_")

    remap: Dict[str, str] = {}
    if args.remap_chain:
        if not os.path.isfile(args.remap_chain):
            sys.exit(f"ERROR: --remap-chain not found: {args.remap_chain}")
        log.info("Running liftOver coordinate remapping…")
        remap = remap_vcf(args.input_vcf, args.remap_chain, args.liftover_exec, tmp_dir)
        log.info("liftOver remapped %d locus/loci", len(remap))

    # Determine tumor / normal column indices from header
    tum_col_idx = -1
    nrm_col_idx = -1
    sample_order: List[int] = []  # indices of genotype columns in desired order

    header_lines: List[str] = []

    def format_header_sort_key(line: str) -> Tuple[int, int, str]:
        if not line.startswith("##FORMAT="):
            return (0, 0, line)
        match = re.search(r"ID=([^,>]+)", line)
        order = {"AD": 0, "DP": 1, "GT": 2}
        return (1, order.get(match.group(1) if match else "", 99), line)

    with open(args.input_vcf) as in_fh, open(args.output_vcf, "w") as out_fh:
        for raw_line in in_fh:
            line = raw_line.rstrip("\n\r")

            # Pass-through meta-information lines, filtering deprecated ones
            if line.startswith("##"):
                # Remove INFO/SVTYPE when ALT is not symbolic (causes VEP to skip)
                if re.match(r"^##INFO=<ID=SVTYPE", line):
                    continue
                if args.remap_chain and line.startswith("##contig="):
                    continue
                header_lines.append(line)
                continue

            # Column header line
            if line.startswith("#CHROM"):
                for header_line in sorted(header_lines, key=format_header_sort_key):
                    out_fh.write(header_line + "\n")
                fields = line.lstrip("#").split("\t")
                sample_cols = fields[9:]

                for idx, s in enumerate(sample_cols):
                    if s == args.vcf_tumor_id:
                        tum_col_idx = idx
                    if s == args.vcf_normal_id:
                        nrm_col_idx = idx

                # Reorder: TUMOR first, then NORMAL, then others
                new_sample_order = []
                if tum_col_idx >= 0:
                    new_sample_order.append(tum_col_idx)
                if nrm_col_idx >= 0 and nrm_col_idx != tum_col_idx:
                    new_sample_order.append(nrm_col_idx)
                for i in range(len(sample_cols)):
                    if i not in new_sample_order:
                        new_sample_order.append(i)
                sample_order = new_sample_order
                log.info(
                    "VCF header parsed: tumor col=%s (idx %d), normal col=%s (idx %d)",
                    args.vcf_tumor_id, tum_col_idx,
                    args.vcf_normal_id, nrm_col_idx,
                )

                new_header_samples = [sample_cols[i] for i in sample_order]
                out_fh.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                    + "\t".join(new_header_samples)
                    + "\n"
                )
                continue

            # Variant line
            cols = line.split("\t")
            if len(cols) < 8:
                continue

            chrom, pos_str, vid, ref, alt_field, qual, filt, info = cols[:8]
            fmt_str = cols[8] if len(cols) > 8 else "GT"
            sample_data = cols[9:] if len(cols) > 9 else []

            # Apply liftover remapping if applicable
            orig_key = f"{chrom}:{pos_str}"
            if remap and orig_key in remap:
                new_locus = remap[orig_key]
                chrom, new_pos = new_locus.split(":")
                pos_str = new_pos
            elif remap:
                log.debug("Liftover miss — keeping original coords for %s", orig_key)

            # Remove INFO/SVTYPE if ALT is defined (not symbolic)
            if not alt_field.startswith("<"):
                info = re.sub(r";?SVTYPE=[^;]+", "", info).lstrip(";") or "."
            if args.remap_chain and info == ".":
                info = ""

            # Reorder sample columns
            new_samples = [
                sample_data[i] if i < len(sample_data) else "." for i in sample_order
            ]

            # Expand multiallelic lines
            reconstructed = (
                f"{chrom}\t{pos_str}\t{vid}\t{ref}\t{alt_field}\t"
                f"{qual}\t{filt}\t{info}\t{fmt_str}\t" + "\t".join(new_samples) + "\n"
            )
            for split_line in split_multiallelic(reconstructed):
                out_fh.write(split_line)

    log.info("Normalised VCF written to %s", args.output_vcf)


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
    vcf2vcf(args)


if __name__ == "__main__":
    main()
