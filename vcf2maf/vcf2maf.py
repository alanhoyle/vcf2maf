#!/usr/bin/env python3
"""
vcf2maf.py  –  Convert a VCF into a MAF, annotating each variant to exactly
one gene isoform.  Mirrors the logic of mskcc/vcf2maf (vcf2maf.pl).

Usage:
    python vcf2maf.py --input-vcf INPUT.vcf --output-maf OUTPUT.maf \\
        --tumor-id TUMOR --normal-id NORMAL

Requires Ensembl VEP to be installed (unless --inhibit-vep is set).
"""

import argparse
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    from .constants import (
        VEP_CONSEQUENCE_PRIORITY,
        VEP_TO_MAF_VARIANT_CLASS,
        MAF_COLUMNS,
        BIOTYPE_PRIORITY,
    )
except ImportError:
    from constants import (  # noqa: E402  (running as a standalone script)
        VEP_CONSEQUENCE_PRIORITY,
        VEP_TO_MAF_VARIANT_CLASS,
        MAF_COLUMNS,
        BIOTYPE_PRIORITY,
    )

log = logging.getLogger("vcf2maf")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Convert a VCF to MAF, annotated via Ensembl VEP.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", default="",
                   help="Config file to load defaults from (default: ~/.vcf2maf.cfg)")
    p.add_argument("--input-vcf", required=True, help="Path to input VCF")
    p.add_argument("--output-maf", required=True, help="Path to output MAF")
    p.add_argument("--tumor-id", default="TUMOR", help="Tumor sample ID")
    p.add_argument("--normal-id", default="NORMAL", help="Normal sample ID")
    p.add_argument(
        "--vcf-tumor-id",
        default=None,
        help="Sample ID in VCF genotype column for tumor (if different from --tumor-id)",
    )
    p.add_argument(
        "--vcf-normal-id",
        default=None,
        help="Sample ID in VCF genotype column for normal (if different from --normal-id)",
    )

    # VEP options
    p.add_argument(
        "--vep-path",
        default=os.path.expanduser("~/miniconda3/bin"),
        help="Directory containing the VEP executable",
    )
    p.add_argument(
        "--vep-data",
        default=os.path.expanduser("~/.vep"),
        help="VEP cache/data directory",
    )
    p.add_argument(
        "--vep-forks", type=int, default=4, help="Number of parallel VEP forks"
    )
    p.add_argument("--buffer-size", type=int, default=5000, help="VEP --buffer_size")
    p.add_argument(
        "--cache-version",
        default="",
        help="VEP offline cache version (e.g. 112); default=auto-detect",
    )
    p.add_argument(
        "--vep-custom", default="", help="Passed to VEP --custom (comma-delimited)"
    )
    p.add_argument("--vep-config", default="", help="Passed to VEP --config")
    p.add_argument(
        "--vep-plugins", default="", help="Passed to VEP --plugin (comma-delimited)"
    )
    p.add_argument(
        "--vep-overwrite",
        action="store_true",
        help="Overwrite pre-existing VEP-annotated VCF",
    )
    p.add_argument(
        "--vep-stats",
        default="",
        help="Path for VEP summary stats file (HTML). "
        "Omit or leave empty to suppress stats output (--no_stats).",
    )

    # Reference / genome
    p.add_argument(
        "--ref-fasta",
        default=os.path.expanduser(
            "~/.vep/homo_sapiens/112_GRCh38/"
            "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        ),
        help="Reference FASTA",
    )
    p.add_argument("--species", default="homo_sapiens")
    p.add_argument(
        "--ncbi-build",
        default="GRCh38",
        help="NCBI reference assembly (GRCh37, GRCh38, GRCm38, …)",
    )

    # MAF options
    p.add_argument(
        "--maf-center", default=".", help="Sequencing center name for MAF Center column"
    )
    p.add_argument(
        "--min-hom-vaf",
        type=float,
        default=0.7,
        help="Minimum VAF to call a variant homozygous",
    )
    p.add_argument(
        "--max-subpop-af",
        type=float,
        default=0.0004,
        help="gnomAD sub-population AF above which variant is tagged common",
    )
    p.add_argument(
        "--retain-info",
        default="",
        help="Comma-separated INFO keys to retain as extra MAF columns",
    )
    p.add_argument(
        "--retain-fmt",
        default="",
        help="Comma-separated FORMAT keys to retain as extra MAF columns",
    )
    p.add_argument(
        "--retain-ann",
        default="",
        help="Comma-separated VEP CSQ fields to retain as extra MAF columns",
    )
    p.add_argument(
        "--custom-enst",
        default="",
        help="File of custom Ensembl transcript IDs (one per line) to prefer",
    )

    # Liftover
    p.add_argument(
        "--remap-chain",
        default="",
        help="UCSC liftOver chain file for coordinate remapping",
    )
    p.add_argument(
        "--liftover-exec", default="liftOver", help="Path to liftOver binary"
    )
    p.add_argument(
        "--samtools",
        default="samtools",
        help="Path to samtools binary (used to fetch flanking reference bases)",
    )

    # Behaviour flags
    p.add_argument(
        "--inhibit-vep",
        action="store_true",
        help="Skip running VEP; parse existing CSQ/ANN annotations if present",
    )
    p.add_argument(
        "--any-allele",
        action="store_true",
        help="Use any ALT allele for annotation, not just tumor ALT",
    )
    p.add_argument(
        "--online",
        action="store_true",
        help="Use VEP REST API instead of offline cache (slow)",
    )
    p.add_argument("--verbose", action="store_true")

    return p


# ---------------------------------------------------------------------------
# VCF parsing helpers
# ---------------------------------------------------------------------------


def open_vcf(path: str):
    """Return a text-mode file handle, handling .gz transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info(info_str: str) -> Dict[str, str]:
    """Parse a VCF INFO field into a dict."""
    d: Dict[str, str] = {}
    for token in info_str.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            d[k] = v
        else:
            d[token] = "1"
    return d


def parse_format(fmt_str: str, sample_str: str) -> Dict[str, str]:
    """Return a dict mapping FORMAT keys to sample values."""
    keys = fmt_str.split(":")
    vals = sample_str.split(":")
    # Pad with '.' if sample has fewer fields than format
    while len(vals) < len(keys):
        vals.append(".")
    return dict(zip(keys, vals))


def get_allele_counts(gt_dict: Dict[str, str]) -> Tuple[int, int, int]:
    """
    Extract (depth, ref_count, alt_count) from a parsed FORMAT dict.
    Handles AD, DP, RD/AD (VarScan), NR/NV, TAR/TIR (Strelka) fields.
    Returns (-1, -1, -1) when counts cannot be determined.
    """
    depth = ref_count = alt_count = -1

    if "DP" in gt_dict:
        try:
            depth = int(gt_dict["DP"])
        except (ValueError, TypeError):
            pass

    # VarScan: RD (ref depth) and AD (alt depth) — must check before standard AD
    # VarScan format has both RD and AD as separate integer fields
    if "RD" in gt_dict and "AD" in gt_dict:
        try:
            ref_count = int(gt_dict["RD"])
            alt_count = int(gt_dict["AD"])
            if depth < 0:
                depth = ref_count + alt_count
        except ValueError:
            pass

    # Standard AD field: REF_count,ALT_count[,ALT2_count,...] (only if VarScan not detected)
    if (
        ref_count < 0
        and "AD" in gt_dict
        and gt_dict["AD"] not in (".", "")
        and "," in gt_dict["AD"]
    ):
        parts = gt_dict["AD"].split(",")
        try:
            ref_count = int(parts[0]) if parts[0].isdigit() else -1
            alt_values = [int(x) for x in parts[1:] if x.isdigit()]
            if alt_values:
                alt_count = sum(alt_values)
            if depth < 0 and ref_count >= 0 and alt_count >= 0:
                depth = ref_count + alt_count
        except (ValueError, IndexError):
            pass

    # NR/NV (Platypus / some others)
    if ref_count < 0 and "NR" in gt_dict and "NV" in gt_dict:
        try:
            total = int(gt_dict["NR"])
            nv = int(gt_dict["NV"])
            alt_count = nv
            ref_count = total - nv
            if depth < 0:
                depth = total
        except ValueError:
            pass

    # Strelka SNVs: xU fields (AU, CU, GU, TU – first value = tier1)
    if alt_count < 0:
        tier1 = {}
        for nuc in "ACGT":
            key = f"{nuc}U"
            if key in gt_dict:
                try:
                    tier1[nuc] = int(gt_dict[key].split(",")[0])
                except (ValueError, IndexError):
                    tier1[nuc] = 0
        if tier1:
            total = sum(tier1.values())
            depth = total
            # alt_count is unknown without knowing the ALT allele here;
            # caller must resolve after knowing ALT base

    return depth, ref_count, alt_count


def determine_variant_type(ref: str, alt: str) -> str:
    """Return MAF Variant_Type: SNP, DNP, TNP, ONP, INS, DEL, or Wildtype."""
    if ref == "-" and alt != "-":
        return "INS"
    if alt == "-" and ref != "-":
        return "DEL"
    if alt in (".", "<DEL>", "<INS>"):
        return "DEL" if alt == "<DEL>" else "INS"
    ref_len = len(ref)
    alt_len = len(alt)
    if ref_len == alt_len:
        diff = sum(1 for a, b in zip(ref, alt) if a != b)
        if diff == 1:
            return "SNP"
        elif diff == 2:
            return "DNP"
        elif diff == 3:
            return "TNP"
        else:
            return "ONP"
    elif alt_len > ref_len:
        return "INS"
    else:
        return "DEL"


def vcf_to_maf_coords(
    chrom: str, pos: int, ref: str, alt: str
) -> Tuple[str, int, int, str, str]:
    """
    Convert VCF-style (1-based, left-aligned) coordinates and alleles to
    MAF-style (1-based inclusive end), normalising indels per MAF spec.

    Returns (chrom, start, end, maf_ref, maf_alt).
    """
    # Strip common prefix (VCF padding for indels)
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    if len(ref) == len(alt):
        # SNP / MNP
        start = pos
        end = pos + len(ref) - 1
        return chrom, start, end, ref, alt

    if len(alt) == 1 and len(ref) > 1:
        # Deletion: VCF uses anchor base; MAF omits it
        start = pos + 1
        end = pos + len(ref) - 1
        maf_ref = ref[1:]
        maf_alt = "-"
        return chrom, start, end, maf_ref, maf_alt

    if len(ref) == 1 and len(alt) > 1:
        # Insertion
        start = pos
        end = pos + 1
        maf_ref = "-"
        maf_alt = alt[1:]
        return chrom, start, end, maf_ref, maf_alt

    # Complex indel – return as-is with MAF coords
    start = pos
    end = pos + max(len(ref), len(alt)) - 1
    return chrom, start, end, ref, alt


# ---------------------------------------------------------------------------
# VEP CSQ annotation parsing
# ---------------------------------------------------------------------------


def parse_csq_header(header_lines: List[str]) -> List[str]:
    """
    Extract the ordered list of CSQ sub-field names from the VCF header line:
      ##INFO=<ID=CSQ,...,Description="...Format: Allele|Consequence|...">
    Returns [] if not found (fall back to ANN if present).
    Matches Perl's [^"']+ pattern to capture field names containing '.' and '/'.
    """
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            m = re.search(r"Format:\s*([^\"']+)[\"']", line)
            if m:
                return [f.strip() for f in m.group(1).split("|")]
    return []


def parse_ann_header(header_lines: List[str]) -> List[str]:
    """
    Extract field order from SnpEff-style ANN field:
      ##INFO=<ID=ANN,...,Description="...Allele|Annotation|...">
    """
    for line in header_lines:
        if line.startswith("##INFO=<ID=ANN"):
            m = re.search(r'"[^"]*\s([\w|]+)"', line)
            if m:
                return m.group(1).split("|")
    return []


def parse_csq_entries(csq_str: str, csq_fields: List[str]) -> List[Dict[str, str]]:
    """
    Parse the CSQ= value from an INFO field into a list of annotation dicts.
    Multiple transcripts are comma-separated; sub-fields are pipe-separated.
    Mirrors Perl: replaces '&' with ',' in every value (s/&/,/g), then applies
    VEP ≥112 field-name aliases (e.g. cDNA.pos/cDNA.length → cDNA_position).
    """
    entries = []
    for entry in csq_str.split(","):
        parts = entry.split("|")
        d: Dict[str, str] = {}
        for i, field in enumerate(csq_fields):
            val = parts[i] if i < len(parts) else ""
            d[field] = val.replace("&", ",")
        # VEP ≥112 renamed several fields; override old names with new values
        for new_name, old_name in [
            ("cDNA.pos/cDNA.length", "cDNA_position"),
            ("CDS.pos/CDS.length", "CDS_position"),
            ("AA.pos/AA.length", "Protein_position"),
            ("HGVS.c", "HGVSc"),
            ("HGVS.p", "HGVSp"),
            ("Distance", "DISTANCE"),
        ]:
            if d.get(new_name):
                d[old_name] = d[new_name]
        entries.append(d)
    return entries


def _transcript_length(csq_entry: Dict[str, str]) -> int:
    """Extract total transcript length from cDNA_position '1480/5971' → 5971."""
    cdna_pos = csq_entry.get("cDNA_position", "")
    if "/" in cdna_pos:
        try:
            return int(cdna_pos.split("/")[-1])
        except (ValueError, IndexError):
            pass
    return 0


def consequence_priority(csq_entry: Dict[str, str]) -> Tuple:
    """
    Return a sort key (consequence_rank, biotype_rank, is_not_canonical,
    -transcript_length) for picking the best transcript annotation.
    Smaller = better; transcript_length is negated so longer sorts first.
    """
    consequence_str = csq_entry.get("Consequence", "")
    # After &→, replacement, consequences are comma-separated
    consequences = consequence_str.split(",")
    rank = min(
        (VEP_CONSEQUENCE_PRIORITY.get(c, 20) for c in consequences),
        default=20,
    )
    biotype = csq_entry.get("BIOTYPE", "")
    bio_rank = BIOTYPE_PRIORITY.get(biotype, 8)
    canonical = 0 if csq_entry.get("CANONICAL", "") == "YES" else 1
    transcript_length = _transcript_length(csq_entry)
    return (rank, bio_rank, canonical, -transcript_length)


def pick_best_csq(
    csq_entries: List[Dict[str, str]],
    custom_enst: Dict[str, bool],
    tumor_allele: str,
    is_deletion: bool = False,
    maf_allele: str = "",
    allele_num: int = 1,
) -> Optional[Dict[str, str]]:
    """
    From a list of per-transcript CSQ annotations, pick the one that best
    describes the effect on the tumor alternate allele, following the same
    priority logic as vcf2maf.pl.
    """
    if not csq_entries:
        return None

    # First restrict to the selected ALT allele when VEP reports ALLELE_NUM.
    # This mirrors upstream vcf2maf.pl and is critical for multiallelic indel
    # records where several alleles may have similar minimized CSQ alleles.
    same_allele = [
        e
        for e in csq_entries
        if not e.get("ALLELE_NUM") or e.get("ALLELE_NUM") == str(allele_num)
    ]
    if not same_allele:
        same_allele = list(csq_entries)

    # Filter to entries matching the tumor allele. With --minimal, VEP may
    # report the minimized MAF-style allele (e.g. A) instead of the anchored
    # VCF ALT (e.g. TA), so try both. VEP can represent deletion alleles as
    # "-", but in multiallelic records a "-" CSQ must not match a non-deletion
    # tumor allele.
    allele_matches = {tumor_allele}
    if maf_allele and maf_allele != "-":
        allele_matches.add(maf_allele)
    matching = [
        e
        for e in same_allele
        if e.get("Allele", "") in allele_matches
        or (is_deletion and e.get("Allele", "") == "-")
    ]

    if not matching:
        matching = list(same_allele)  # fall back to selected allele entries

    # Prefer custom ENST transcripts first
    if custom_enst:
        custom_matching = [e for e in matching if e.get("Feature", "") in custom_enst]
        if custom_matching:
            matching = custom_matching

    # Sort by priority; ties broken by ENST ID for reproducibility
    matching.sort(key=lambda e: consequence_priority(e))
    return matching[0]


# ---------------------------------------------------------------------------
# Variant classification
# ---------------------------------------------------------------------------


def classify_variant(csq: Dict[str, str], variant_type: str, ref: str, alt: str) -> str:
    """
    Map VEP consequence(s) to a single MAF Variant_Classification string.
    """
    consequences = csq.get("Consequence", "").split(",")
    if variant_type == "INS" and "protein_altering_variant" in consequences:
        return "In_Frame_Ins"
    if variant_type == "DEL" and "protein_altering_variant" in consequences:
        return "In_Frame_Del"
    # Use the highest-priority (lowest rank) consequence
    best_cons = min(
        consequences,
        key=lambda c: VEP_CONSEQUENCE_PRIORITY.get(c, 20),
    )
    maf_class = VEP_TO_MAF_VARIANT_CLASS.get(best_cons, "")

    # Refine Frame_Shift based on variant type
    if maf_class == "Frame_Shift":
        maf_class = "Frame_Shift_Ins" if variant_type == "INS" else "Frame_Shift_Del"

    # Refine Silent for splice region
    if best_cons == "splice_region_variant":
        exon = csq.get("EXON", "")
        intron = csq.get("INTRON", "")
        if exon:
            maf_class = "Splice_Region"
        elif intron:
            maf_class = "Splice_Region"

    if not maf_class:
        # Fall back based on position
        if variant_type == "INS":
            maf_class = "In_Frame_Ins"
        elif variant_type == "DEL":
            maf_class = "In_Frame_Del"
        else:
            maf_class = "Missense_Mutation"

    return maf_class


def hgvsp_short(hgvsp_long: str) -> str:
    """
    Convert HGVSp (3-letter amino-acid codes) to short 1-letter form.
    e.g. p.Glu746_Ala750del  ->  p.E746_A750del
         p.Thr125=           ->  p.T125=
    """
    if not hgvsp_long or "%" in hgvsp_long:
        # URL-encode artifacts – just return as-is
        return hgvsp_long

    aa3to1 = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Ser": "S",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
        "Ter": "*",
        "Sec": "U",
        "Pyl": "O",
        "Xaa": "X",
        "Xle": "J",
    }

    def replace_aa(match: re.Match) -> str:
        return aa3to1.get(match.group(0), match.group(0))

    short = re.sub(r"[A-Z][a-z]{2}", replace_aa, hgvsp_long)
    return short


# ---------------------------------------------------------------------------
# VEP runner
# ---------------------------------------------------------------------------


def run_vep(input_vcf: str, vep_vcf: str, args: argparse.Namespace) -> None:
    """Build and execute the VEP command."""
    vep_bin = os.path.join(args.vep_path, "vep")
    if not os.path.isfile(vep_bin):
        vep_bin = shutil.which("vep") or "vep"

    cmd = [
        vep_bin,
        "--input_file",
        input_vcf,
        "--output_file",
        vep_vcf,
        "--format",
        "vcf",
        "--vcf",
        "--everything",
        "--allele_number",
        # Extra flags used by vcf2maf.pl but not covered by --everything:
        "--total_length",  # cDNA_position as pos/total (e.g. 1480/5971)
        "--xref_refseq",  # adds RefSeq column to CSQ
        "--flag_pick_allele",  # adds PICK flag to chosen transcript in CSQ
        "--pick_order",
        "canonical,tsl,biotype,rank,ccds,length",
        "--check_existing",
        "--failed",
        "1",
        "--no_escape",
        "--shift_hgvs",
        "1",
        "--cache",
        "--offline",
        "--dir_cache",
        args.vep_data,
        "--fasta",
        args.ref_fasta,
        "--minimal",
        "--fork",
        str(args.vep_forks),
        "--buffer_size",
        str(args.buffer_size),
        "--species",
        args.species,
        "--assembly",
        args.ncbi_build,
    ]
    # Stats file: write to the requested path, or suppress entirely
    if getattr(args, "vep_stats", ""):
        cmd += ["--stats_file", args.vep_stats]
    else:
        cmd += ["--no_stats"]
    if args.cache_version:
        cmd += ["--cache_version", args.cache_version]
    if args.vep_custom:
        cmd += ["--custom", args.vep_custom]
    if args.vep_config:
        cmd += ["--config", args.vep_config]
    if args.vep_plugins:
        for plugin in args.vep_plugins.split(","):
            cmd += ["--plugin", plugin.strip()]
    if getattr(args, "vep_overwrite", False):
        cmd += ["--force_overwrite"]

    log.info("Running VEP: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=not args.verbose)
    if result.returncode != 0:
        stderr = result.stderr.decode() if result.stderr else ""
        raise RuntimeError(f"VEP failed (exit {result.returncode}):\n{stderr}")


# ---------------------------------------------------------------------------
# Main conversion logic
# ---------------------------------------------------------------------------


def vcf2maf(args: argparse.Namespace) -> None:
    """Full VCF → MAF conversion pipeline."""

    input_vcf = args.input_vcf
    output_maf = args.output_maf

    # Basic validation
    if not os.path.isfile(input_vcf):
        sys.exit(f"ERROR: --input-vcf not found: {input_vcf}")
    if not os.path.isfile(args.ref_fasta):
        sys.exit(f"ERROR: --ref-fasta not found: {args.ref_fasta}")
    if input_vcf.endswith((".gz", ".bz2", ".bcf")):
        sys.exit("ERROR: --input-vcf cannot be compressed. Please decompress first.")

    # Load Ensembl → Entrez gene ID map (mirrors Perl $script_dir/data/… lookup)
    entrez_map = _load_entrez_map()

    # Pre-fetch flanking reference bases for all variants (mirrors Perl pre-VEP step)
    samtools = getattr(args, "samtools", "samtools")
    flanking_map = _fetch_all_flanking_bps(input_vcf, args.ref_fasta, samtools)

    # Load custom transcript list
    custom_enst: Dict[str, bool] = {}
    if args.custom_enst:
        if not os.path.isfile(args.custom_enst):
            sys.exit(f"ERROR: --custom-enst file not found: {args.custom_enst}")
        with open(args.custom_enst) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    custom_enst[line.split()[0]] = True
        log.info("Loaded %d custom ENST IDs", len(custom_enst))

    # Determine which column in VCF corresponds to tumor/normal
    vcf_tumor_id = args.vcf_tumor_id or args.tumor_id
    vcf_normal_id = args.vcf_normal_id or args.normal_id

    # -----------------------------------------------------------------------
    # Step 1: Run VEP (or reuse existing .vep.vcf)
    # -----------------------------------------------------------------------
    vep_vcf = re.sub(r"\.vcf$", ".vep.vcf", input_vcf)
    if vep_vcf == input_vcf:
        vep_vcf = input_vcf + ".vep.vcf"

    if not args.inhibit_vep:
        if os.path.isfile(vep_vcf) and not args.vep_overwrite:
            log.info(
                "Found existing VEP-annotated VCF: %s (use --vep-overwrite to rerun)",
                vep_vcf,
            )
        else:
            run_vep(input_vcf, vep_vcf, args)
    else:
        vep_vcf = input_vcf  # parse annotations from the input itself

    # -----------------------------------------------------------------------
    # Step 2: Parse VEP-annotated VCF
    # -----------------------------------------------------------------------
    retain_info_keys = [k for k in args.retain_info.split(",") if k]
    retain_fmt_keys = [k for k in args.retain_fmt.split(",") if k]
    retain_ann_keys = [k for k in args.retain_ann.split(",") if k]

    extra_cols = retain_info_keys + retain_fmt_keys + retain_ann_keys

    header_lines: List[str] = []
    csq_fields: List[str] = []
    ann_fields: List[str] = []
    tum_col_idx = -1
    nrm_col_idx = -1
    sample_cols: List[str] = []

    maf_rows: List[Dict] = []

    with open_vcf(vep_vcf) as vcf_fh:
        for raw_line in vcf_fh:
            line = raw_line.rstrip("\n\r")

            # ---- Header -----------------------------------------------
            if line.startswith("##"):
                header_lines.append(line)
                continue

            if line.startswith("#CHROM"):
                header_lines.append(line)
                fields = line.lstrip("#").split("\t")
                # Columns 9+ are sample genotype columns
                sample_cols = fields[9:]
                # Locate tumor / normal columns
                for idx, s in enumerate(sample_cols):
                    if s == vcf_tumor_id:
                        tum_col_idx = idx
                    if s == vcf_normal_id:
                        nrm_col_idx = idx

                csq_fields = parse_csq_header(header_lines)
                if not csq_fields:
                    ann_fields = parse_ann_header(header_lines)
                continue

            # ---- Variant line -----------------------------------------
            cols = line.split("\t")
            if len(cols) < 8:
                continue

            chrom, pos_str, vcf_id, ref, alt_field, qual, filter_str, info_str = cols[
                :8
            ]
            pos = int(pos_str)

            fmt_str = cols[8] if len(cols) > 8 else ""
            sample_data = cols[9:] if len(cols) > 9 else []

            info = parse_info(info_str)

            # Parse genotype dicts for tumor and normal
            tum_gt = {}
            nrm_gt = {}
            if fmt_str and tum_col_idx >= 0 and tum_col_idx < len(sample_data):
                tum_gt = parse_format(fmt_str, sample_data[tum_col_idx])
            if fmt_str and nrm_col_idx >= 0 and nrm_col_idx < len(sample_data):
                nrm_gt = parse_format(fmt_str, sample_data[nrm_col_idx])

            # Determine tumor ALT allele from genotype
            alts = alt_field.split(",")
            tumor_allele = _resolve_tumor_allele(tum_gt, ref, alts)

            # Track 1-based allele index (mirrors Perl $var_allele_idx into @alleles)
            var_allele_idx = next(
                (i for i, a in enumerate(alts, 1) if a == tumor_allele), 1
            )

            # Skip if tumor allele is the reference
            if tumor_allele == ref and not args.any_allele:
                continue

            # Convert to MAF coordinates
            maf_chrom, start, end, maf_ref, maf_alt = vcf_to_maf_coords(
                chrom, pos, ref, tumor_allele
            )
            variant_type = determine_variant_type(maf_ref, maf_alt)

            # ---- Parse CSQ / ANN annotations --------------------------
            csq_entries: List[Dict[str, str]] = []
            if "CSQ" in info and csq_fields:
                csq_entries = parse_csq_entries(info["CSQ"], csq_fields)
            elif "ANN" in info and ann_fields:
                csq_entries = parse_csq_entries(info["ANN"], ann_fields)

            is_deletion_allele = maf_alt == "-"
            best_csq = pick_best_csq(
                csq_entries,
                custom_enst,
                tumor_allele,
                is_deletion_allele,
                maf_alt,
                var_allele_idx,
            )
            if best_csq is None:
                best_csq = {}

            # ---- Variant classification --------------------------------
            var_class = classify_variant(best_csq, variant_type, maf_ref, maf_alt)

            # ---- Allele counts ----------------------------------------
            t_depth, t_ref_count, t_alt_count = get_allele_counts(tum_gt)
            n_depth, n_ref_count, n_alt_count = get_allele_counts(nrm_gt)
            # Raw AD string values — preserves "." per Perl @tum_depths logic
            t_ref_str, t_alt_str = _ad_str_counts(tum_gt, var_allele_idx)
            n_ref_str, n_alt_str = _ad_str_counts(nrm_gt, var_allele_idx)

            # ---- Zygosity / allele assignment --------------------------
            tum_seq_allele1, tum_seq_allele2 = _assign_alleles(
                tum_gt, ref, maf_ref, maf_alt, t_alt_count, t_depth, args.min_hom_vaf
            )
            # Build full allele list so GT indices > 1 map to the right MAF allele
            full_allele_list = [maf_ref]
            for _vcf_alt in alts:
                _, _, _, _, _extra_maf_alt = vcf_to_maf_coords(chrom, pos, ref, _vcf_alt)
                full_allele_list.append(_extra_maf_alt)
            nrm_seq_allele1, nrm_seq_allele2 = _assign_normal_alleles(
                nrm_gt, maf_ref, maf_alt, full_allele_list
            )

            # ---- FILTER / common variant tagging ----------------------
            common_tag = _is_common_variant(best_csq, args.max_subpop_af)
            # Mirrors Perl: keep VCF FILTER as-is; only default empty to "."
            filter_tag = filter_str if filter_str else "."
            if common_tag and "common_variant" not in filter_tag:
                filter_tag = (
                    "common_variant"
                    if filter_tag in (".", "PASS")
                    else filter_tag + ";common_variant"
                )

            # ---- HGVSc / HGVSp prefix stripping (mirrors Perl lines 787-791) --
            # VEP ≥112 prefixes with transcript/protein ID: "ENST…:c.xxx", "ENSP…:p.xxx"
            hgvsc_raw = best_csq.get("HGVSc", "")
            hgvsp_raw = best_csq.get("HGVSp", "")
            hgvsc_clean = re.sub(r"^.*:", "", hgvsc_raw) if hgvsc_raw else ""
            hgvsp_clean = re.sub(r"^.*:", "", hgvsp_raw) if hgvsp_raw else ""
            # Also strip "c.xxx(p.xxx)" wrapping
            if hgvsp_clean and hgvsp_clean.startswith("c."):
                m = re.search(r"\((p\.\S+)\)", hgvsp_clean)
                if m:
                    hgvsp_clean = m.group(1)
            hgvsp_s = hgvsp_short(hgvsp_clean)

            # Synthetic HGVSp_Short for splice acceptor/donor variants (mirrors Perl lines 824-830)
            best_consequence = best_csq.get("Consequence", "").split(",")[0]
            if best_consequence in ("splice_acceptor_variant", "splice_donor_variant"):
                m_cpos = re.match(r"^c\.(\d+)", hgvsc_clean)
                if m_cpos:
                    c_pos = max(1, int(m_cpos.group(1)))
                    p_pos = round((c_pos + c_pos % 3) / 3)
                    hgvsp_s = f"p.X{p_pos}_splice"

            # ---- Collect all transcript effects for all_effects column -
            all_effects = _build_all_effects(
                csq_entries, tumor_allele, is_deletion_allele, maf_alt, var_allele_idx
            )

            # ---- Build MAF row ----------------------------------------
            ensg_id = best_csq.get("Gene", "")
            row: Dict = {
                "Hugo_Symbol": best_csq.get("SYMBOL", "Unknown"),
                "Entrez_Gene_Id": entrez_map.get(ensg_id, "0"),
                "Center": args.maf_center,
                "NCBI_Build": args.ncbi_build,
                "Chromosome": maf_chrom,
                "Start_Position": str(start),
                "End_Position": str(end),
                "Strand": "+",
                "Variant_Classification": var_class,
                "Variant_Type": variant_type,
                "Reference_Allele": maf_ref,
                "Tumor_Seq_Allele1": tum_seq_allele1,
                "Tumor_Seq_Allele2": tum_seq_allele2,
                "dbSNP_RS": _extract_dbsnp(best_csq.get("Existing_variation", "")),
                "dbSNP_Val_Status": "",
                "Tumor_Sample_Barcode": args.tumor_id,
                "Matched_Norm_Sample_Barcode": args.normal_id,
                "Match_Norm_Seq_Allele1": nrm_seq_allele1,
                "Match_Norm_Seq_Allele2": nrm_seq_allele2,
                "Tumor_Validation_Allele1": "",
                "Tumor_Validation_Allele2": "",
                "Match_Norm_Validation_Allele1": "",
                "Match_Norm_Validation_Allele2": "",
                "Verification_Status": "",
                "Validation_Status": "",
                "Mutation_Status": "",
                "Sequencing_Phase": "",
                "Sequence_Source": "",
                "Validation_Method": "",
                "Score": "",
                "BAM_File": "",
                "Sequencer": "",
                "Tumor_Sample_UUID": "",
                "Matched_Norm_Sample_UUID": "",
                "HGVSc": hgvsc_clean,
                "HGVSp": hgvsp_clean,
                "HGVSp_Short": hgvsp_s,
                "Transcript_ID": best_csq.get("Feature", ""),
                "Exon_Number": best_csq.get("EXON", ""),
                "t_depth": str(t_depth) if t_depth >= 0 else "",
                "t_ref_count": (
                    t_ref_str
                    if t_ref_str != ""
                    else (str(t_ref_count) if t_ref_count >= 0 else "")
                ),
                "t_alt_count": (
                    t_alt_str
                    if t_alt_str != ""
                    else (str(t_alt_count) if t_alt_count >= 0 else "")
                ),
                "n_depth": str(n_depth) if n_depth >= 0 else "",
                "n_ref_count": (
                    n_ref_str
                    if n_ref_str != ""
                    else (str(n_ref_count) if n_ref_count >= 0 else "")
                ),
                "n_alt_count": (
                    n_alt_str
                    if n_alt_str != ""
                    else (str(n_alt_count) if n_alt_count >= 0 else "")
                ),
                "all_effects": all_effects,
                # VEP pass-through columns
                "Allele": best_csq.get("Allele", ""),
                "Gene": best_csq.get("Gene", ""),
                "Feature": best_csq.get("Feature", ""),
                "Feature_type": best_csq.get("Feature_type", ""),
                "One_Consequence": best_csq.get("Consequence", "").split("&")[0],
                "Consequence": best_csq.get("Consequence", ""),
                "cDNA_position": best_csq.get("cDNA_position", ""),
                "CDS_position": best_csq.get("CDS_position", ""),
                "Protein_position": best_csq.get("Protein_position", ""),
                "Amino_acids": best_csq.get("Amino_acids", ""),
                "Codons": best_csq.get("Codons", ""),
                "Existing_variation": best_csq.get("Existing_variation", ""),
                "ALLELE_NUM": best_csq.get("ALLELE_NUM", ""),
                "DISTANCE": best_csq.get("DISTANCE", ""),
                "TRANSCRIPT_STRAND": best_csq.get("STRAND", ""),
                "SYMBOL": best_csq.get("SYMBOL", ""),
                "SYMBOL_SOURCE": best_csq.get("SYMBOL_SOURCE", ""),
                "HGNC_ID": best_csq.get("HGNC_ID", ""),
                "BIOTYPE": best_csq.get("BIOTYPE", ""),
                "CANONICAL": best_csq.get("CANONICAL", ""),
                "CCDS": best_csq.get("CCDS", ""),
                "ENSP": best_csq.get("ENSP", ""),
                "SWISSPROT": best_csq.get("SWISSPROT", ""),
                "TREMBL": _strip_organism_suffix(best_csq.get("TREMBL", "")),
                "UNIPARC": best_csq.get("UNIPARC", ""),
                "RefSeq": best_csq.get("RefSeq", ""),
                "SIFT": best_csq.get("SIFT", ""),
                "PolyPhen": best_csq.get("PolyPhen", ""),
                "EXON": best_csq.get("EXON", ""),
                "INTRON": best_csq.get("INTRON", ""),
                "DOMAINS": best_csq.get("DOMAINS", ""),
                "gnomAD_AF": best_csq.get("gnomAD_AF", "")
                or best_csq.get("gnomADe_AF", ""),
                "gnomAD_AFR_AF": best_csq.get("gnomAD_AFR_AF", "")
                or best_csq.get("gnomADe_AFR_AF", ""),
                "gnomAD_AMR_AF": best_csq.get("gnomAD_AMR_AF", "")
                or best_csq.get("gnomADe_AMR_AF", ""),
                "gnomAD_ASJ_AF": best_csq.get("gnomAD_ASJ_AF", "")
                or best_csq.get("gnomADe_ASJ_AF", ""),
                "gnomAD_EAS_AF": best_csq.get("gnomAD_EAS_AF", "")
                or best_csq.get("gnomADe_EAS_AF", ""),
                "gnomAD_FIN_AF": best_csq.get("gnomAD_FIN_AF", "")
                or best_csq.get("gnomADe_FIN_AF", ""),
                "gnomAD_NFE_AF": best_csq.get("gnomAD_NFE_AF", "")
                or best_csq.get("gnomADe_NFE_AF", ""),
                "gnomAD_OTH_AF": best_csq.get("gnomAD_OTH_AF", "")
                or best_csq.get("gnomADe_OTH_AF", ""),
                "gnomAD_SAS_AF": best_csq.get("gnomAD_SAS_AF", "")
                or best_csq.get("gnomADe_SAS_AF", ""),
                "MAX_AF": best_csq.get("MAX_AF", ""),
                "MAX_AF_POPS": best_csq.get("MAX_AF_POPS", ""),
                "gnomADe_AF": best_csq.get("gnomADe_AF", ""),
                "gnomADe_AFR_AF": best_csq.get("gnomADe_AFR_AF", ""),
                "gnomADe_AMR_AF": best_csq.get("gnomADe_AMR_AF", ""),
                "gnomADe_ASJ_AF": best_csq.get("gnomADe_ASJ_AF", ""),
                "gnomADe_EAS_AF": best_csq.get("gnomADe_EAS_AF", ""),
                "gnomADe_FIN_AF": best_csq.get("gnomADe_FIN_AF", ""),
                "gnomADe_NFE_AF": best_csq.get("gnomADe_NFE_AF", ""),
                "gnomADe_OTH_AF": best_csq.get("gnomADe_OTH_AF", ""),
                "gnomADe_SAS_AF": best_csq.get("gnomADe_SAS_AF", ""),
                "CLIN_SIG": best_csq.get("CLIN_SIG", ""),
                "SOMATIC": best_csq.get("SOMATIC", ""),
                "PUBMED": best_csq.get("PUBMED", ""),
                "TRANSCRIPTION_FACTORS": best_csq.get("TRANSCRIPTION_FACTORS", ""),
                "MOTIF_NAME": best_csq.get("MOTIF_NAME", ""),
                "MOTIF_POS": best_csq.get("MOTIF_POS", ""),
                "HIGH_INF_POS": best_csq.get("HIGH_INF_POS", ""),
                "MOTIF_SCORE_CHANGE": best_csq.get("MOTIF_SCORE_CHANGE", ""),
                "IMPACT": best_csq.get("IMPACT", ""),
                "PICK": best_csq.get("PICK", ""),
                "VARIANT_CLASS": best_csq.get("VARIANT_CLASS", ""),
                "TSL": best_csq.get("TSL", ""),
                "HGVS_OFFSET": best_csq.get("HGVS_OFFSET", ""),
                "PHENO": best_csq.get("PHENO", ""),
                "GENE_PHENO": best_csq.get("GENE_PHENO", ""),
                "FILTER": filter_tag,
                "flanking_bps": flanking_map.get((chrom, pos), ""),
                "variant_id": vcf_id,
                "variant_qual": qual,
                "vcf_id": vcf_id,
                "vcf_qual": qual,
            }

            # Retain extra INFO fields
            for k in retain_info_keys:
                row[k] = info.get(k, "")

            # Retain extra FORMAT fields (tumor sample)
            for k in retain_fmt_keys:
                row[k] = tum_gt.get(k, "")

            # Retain extra VEP annotation fields
            for k in retain_ann_keys:
                row[k] = best_csq.get(k, "")

            maf_rows.append(row)

    # -----------------------------------------------------------------------
    # Step 3: Write MAF
    # -----------------------------------------------------------------------
    output_cols = list(MAF_COLUMNS) + extra_cols
    # Deduplicate while preserving order
    seen: set = set()
    final_cols = []
    for c in output_cols:
        if c not in seen:
            seen.add(c)
            final_cols.append(c)

    with open(output_maf, "w") as maf_fh:
        maf_fh.write("#version 2.4\n")
        maf_fh.write("\t".join(final_cols) + "\n")
        for row in maf_rows:
            maf_fh.write("\t".join(str(row.get(c, "")) for c in final_cols) + "\n")

    log.info("Wrote %d variants to %s", len(maf_rows), output_maf)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _resolve_tumor_allele(tum_gt: Dict[str, str], ref: str, alts: List[str]) -> str:
    """
    Given a parsed tumor FORMAT dict, return the tumor ALT allele.
    Falls back to alts[0] if the genotype is absent or ambiguous.
    """
    gt_val = tum_gt.get("GT", "")
    if gt_val in (".", "./.", ".|."):
        return alts[0] if alts else ref

    # Parse diploid genotype: 0/1, 1/1, 0/2, etc.
    allele_indices = re.split(r"[/|]", gt_val)
    non_ref = [i for i in allele_indices if i not in ("0", ".", "")]
    if non_ref:
        idx = int(non_ref[0]) - 1  # 1-based → 0-based into alts
        if 0 <= idx < len(alts):
            return alts[idx]

    return alts[0] if alts else ref


def _assign_alleles(
    gt_dict: Dict[str, str],
    ref: str,
    maf_ref: str,
    maf_alt: str,
    alt_count: int,
    depth: int,
    min_hom_vaf: float,
) -> Tuple[str, str]:
    """
    Tumor allele assignment: Allele2 is always the variant; Allele1 is the
    first non-variant allele from GT (mirrors Perl vcf2maf.pl logic).
    """
    gt_val = gt_dict.get("GT", "")
    allele_indices = (
        re.split(r"[/|]", gt_val) if gt_val and gt_val not in (".", "./.") else []
    )

    if maf_ref == "-":
        return maf_ref, maf_alt

    # Homozygous alt?
    if set(allele_indices) == {"1"} or (
        depth > 0 and alt_count >= 0 and alt_count / depth >= min_hom_vaf
    ):
        return maf_alt, maf_alt

    # Heterozygous or unknown: allele2 must be the variant
    return maf_ref, maf_alt


def _assign_normal_alleles(
    gt_dict: Dict[str, str], maf_ref: str, maf_alt: str,
    allele_list: Optional[List[str]] = None,
) -> Tuple[str, str]:
    """
    Normal allele assignment: read both GT indices directly and map to
    ref/alt (0=ref, 1=alt).  Defaults to ref/ref when GT is absent.
    Mirrors Perl: defaults to $ref/$ref then overrides from GT indices.

    allele_list: index 0 = maf_ref, index 1 = alts[0] in MAF coords,
    index 2 = alts[1] in MAF coords, etc.  When None, falls back to
    [maf_ref, maf_alt] (the simple biallelic case).
    """
    gt_val = gt_dict.get("GT", "")
    if not gt_val or gt_val in (".", "./."):
        return maf_ref, maf_ref

    indices = re.split(r"[/|]", gt_val)
    if allele_list is None:
        allele_list = [maf_ref, maf_alt]

    def idx_to_allele(s: str) -> str:
        try:
            i = int(s)
            return allele_list[i] if i < len(allele_list) else maf_ref
        except (ValueError, TypeError):
            return maf_ref

    a1 = idx_to_allele(indices[0])
    a2 = idx_to_allele(indices[-1])
    return a1, a2


def _extract_dbsnp(existing_variation: str) -> str:
    """Extract the first rs ID from the Existing_variation field, or 'novel'."""
    for token in re.split(r"[&,]", existing_variation) if existing_variation else []:
        if re.match(r"^rs\d+$", token):
            return token
    return "novel"


def _is_common_variant(csq: Dict[str, str], max_af: float) -> bool:
    """Return True if any gnomAD sub-population AF exceeds max_af."""
    af_keys = [
        "gnomAD_AF",
        "gnomAD_AFR_AF",
        "gnomAD_AMR_AF",
        "gnomAD_ASJ_AF",
        "gnomAD_EAS_AF",
        "gnomAD_FIN_AF",
        "gnomAD_NFE_AF",
        "gnomAD_OTH_AF",
        "gnomAD_SAS_AF",
        "gnomADe_AF",
        "AF",
    ]
    for k in af_keys:
        val = csq.get(k, "")
        if val and val != ".":
            try:
                if float(val) > max_af:
                    return True
            except ValueError:
                pass
    return False


def _load_entrez_map() -> Dict[str, str]:
    """
    Load Ensembl-gene-ID → Entrez-ID map from the bundled data file.
    Mirrors Perl: $script_dir/data/ensg_to_entrez_id_map_ensembl_feb2014.tsv
    """
    candidates = [
        Path(__file__).parent.parent
        / "data"
        / "ensg_to_entrez_id_map_ensembl_feb2014.tsv",
        Path(__file__).parent / "data" / "ensg_to_entrez_id_map_ensembl_feb2014.tsv",
        Path("/opt/data/ensg_to_entrez_id_map_ensembl_feb2014.tsv"),
    ]
    for p in candidates:
        if p.exists():
            mapping: Dict[str, str] = {}
            with open(p) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2:
                        mapping[parts[0]] = parts[1]
            return mapping
    return {}


def _strip_organism_suffix(value: str) -> str:
    """
    Strip trailing _ORGANISM suffix (e.g. _HUMAN) that VEP ≥112 adds to
    TrEMBL accessions (Q96QW8_HUMAN → Q96QW8).  Comma-separated entries
    are each stripped individually.
    """
    if not value:
        return value
    parts = [re.sub(r"_[A-Z0-9]+$", "", p) for p in value.split(",") if p]
    return ",".join(parts)


def _ad_str_counts(gt_dict: Dict[str, str], var_allele_idx: int) -> Tuple[str, str]:
    """
    Return (ref_count_str, alt_count_str) by reading the AD field directly,
    preserving '.' values as-is.  Mirrors Perl: @tum_depths = split(",", AD).
    """
    # VarScan: separate RD (ref) and AD (alt) integer fields
    if "RD" in gt_dict and "AD" in gt_dict and "," not in gt_dict.get("AD", ""):
        return gt_dict.get("RD", ""), gt_dict.get("AD", "")
    ad = gt_dict.get("AD", "")
    if not ad:
        return "", ""
    parts = ad.split(",")
    ref_str = parts[0] if parts else ""
    alt_str = parts[var_allele_idx] if var_allele_idx < len(parts) else ""
    return ref_str, alt_str


def _fetch_all_flanking_bps(
    input_vcf: str, ref_fasta: str, samtools: str = "samtools"
) -> Dict[Tuple[str, int], str]:
    """
    Pre-fetch flanking reference bases for every variant in input_vcf.
    For each variant at (chrom, pos) with ref allele of length L, fetches
    the region chrom:(pos-1)-(pos+L) — 1 base before + ref + 1 base after.
    Returns dict keyed by (chrom, pos) → uppercase base string.
    Mirrors Perl's pre-VEP samtools faidx batch fetch.
    """
    region_to_key: Dict[str, Tuple[str, int]] = {}
    try:
        with open_vcf(input_vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 4:
                    continue
                chrom, pos_str, _id, ref = cols[0], cols[1], cols[2], cols[3]
                try:
                    pos = int(pos_str)
                except ValueError:
                    continue
                region = f"{chrom}:{pos - 1}-{pos + len(ref)}"
                region_to_key[region] = (chrom, pos)
    except OSError:
        return {}

    if not region_to_key:
        return {}

    cmd = [samtools, "faidx", ref_fasta] + list(region_to_key)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    except (OSError, FileNotFoundError):
        return {}
    if result.returncode != 0:
        return {}

    flanking: Dict[Tuple[str, int], str] = {}
    current_region: Optional[str] = None
    seq_parts: List[str] = []

    def _store(region: Optional[str], parts: List[str]) -> None:
        if region and parts:
            key = region_to_key.get(region)
            if key:
                flanking[key] = "".join(parts).upper()

    for line in result.stdout.splitlines():
        if line.startswith(">"):
            _store(current_region, seq_parts)
            current_region = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    _store(current_region, seq_parts)

    return flanking


def _build_all_effects(
    csq_entries: List[Dict[str, str]],
    tumor_allele: str,
    is_deletion: bool = False,
    maf_allele: str = "",
    allele_num: int = 1,
) -> str:
    """
    Build the all_effects column.  Mirrors Perl vcf2maf.pl:
    - Sort by biotype priority, consequence priority, transcript_length desc
    - Format: 'Gene,Consequence,HGVSp,Transcript_ID,RefSeq;' per entry
    """
    allele_matches = {tumor_allele}
    if maf_allele and maf_allele != "-":
        allele_matches.add(maf_allele)
    same_allele = [
        e
        for e in csq_entries
        if not e.get("ALLELE_NUM") or e.get("ALLELE_NUM") == str(allele_num)
    ]
    if not same_allele:
        same_allele = list(csq_entries)
    matching = [
        e
        for e in same_allele
        if (
            e.get("Allele", "") in allele_matches
            or e.get("Allele", "") == ""
            or (is_deletion and e.get("Allele", "") == "-")
        )
    ]
    if not matching:
        matching = list(same_allele)

    def _sort_key(e: Dict[str, str]) -> tuple:
        cons_str = e.get("Consequence", "")
        cons_list = [c for c in cons_str.split(",") if c]
        rank = min((VEP_CONSEQUENCE_PRIORITY.get(c, 20) for c in cons_list), default=20)
        bio_rank = BIOTYPE_PRIORITY.get(e.get("BIOTYPE", ""), 8)
        return (bio_rank, rank, -_transcript_length(e))

    matching.sort(key=_sort_key)

    out = ""
    for e in matching:
        gene = e.get("SYMBOL", "")
        cons_str = e.get("Consequence", "")
        cons = cons_str.split(",")[0] if cons_str else ""
        hgvsp = re.sub(r"^.*:", "", e.get("HGVSp", ""))
        if hgvsp.startswith("c."):
            m_paren = re.search(r"\((p\.\S+)\)", hgvsp)
            if m_paren:
                hgvsp = m_paren.group(1)
        tid = e.get("Feature", "")
        refseq = e.get("RefSeq", "")
        if cons and tid:
            out += f"{gene},{cons},{hgvsp},{tid},{refseq};"
    return out


# ---------------------------------------------------------------------------
# CLI entry point
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
    vcf2maf(args)


if __name__ == "__main__":
    main()
