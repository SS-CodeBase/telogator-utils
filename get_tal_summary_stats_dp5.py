#!/usr/bin/env python3

"""
Low-Stringency Telomere Length Summary (Depth >= 5)
---------------------------------------------------
This variant of the summarization script applies a low-stringency 
quality control threshold. It enforces a minimum support requirement 
(>= 5 reads) and ranks alleles based on a combination of read depth, 
clustering status, and maximum p75 length.
"""

import pandas as pd
import numpy as np
import argparse
import os
import sys
import re

def natural_sort(values):
    """Sorts alphanumeric strings naturally (e.g., chr2 before chr10)."""
    def key(s):
        return [int(x) if x.isdigit() else x.lower()
                for x in re.split(r'(\d+)', str(s))]
    return sorted(values, key=key)

def max_read_tl(tl_string):
    """Parses a comma-separated string of read lengths and returns the max value."""
    try:
        vals = [int(x) for x in str(tl_string).split(",") if x.strip()]
        return max(vals) if vals else np.nan
    except Exception:
        return np.nan

def load_chr_arms(species_table, species):
    """Loads chromosome mapping and generates expected p and q arms."""
    df = pd.read_csv(species_table, sep="\t")
    df = df[df["REF_ACC"] == species]

    if df.empty:
        sys.exit(f"Species not found in chromosome table: {species}")

    arms = []
    arm_to_gb = {}
    arm_to_rs = {}

    for _, r in df.iterrows():
        chrom = r["CHROM"]
        for suffix in ["p", "q"]:
            arm = f"{chrom}{suffix}"
            arms.append(arm)
            arm_to_gb[arm] = r["GB_ACC"]
            arm_to_rs[arm] = r["RS_ACC"]

    return natural_sort(arms), arm_to_gb, arm_to_rs

def process_telomere_data(input_file, species_table, species, output_file, error_file, allele_mode="collapsed"):
    arms, arm_to_gb, arm_to_rs = load_chr_arms(species_table, species)

    # 1. Handle missing inputs gracefully
    if not os.path.exists(input_file):
        with open(error_file, "w") as f:
            f.write(f"{species}: tlens file not found\n")
        
        base = pd.DataFrame({"chr_arm": arms})
        base["chr_gb_acc"] = base["chr_arm"].map(arm_to_gb)
        base["chr_rs_acc"] = base["chr_arm"].map(arm_to_rs)
        
        if allele_mode == "haplotype":
            base = base.loc[base.index.repeat(2)]
            base["haplotype"] = [1, 2] * (len(base) // 2)
            
        base["telomere_length_based_on_TL_p75"] = np.nan
        base.to_csv(output_file, sep="\t", index=False, na_rep="n/a")
        return

    # 2. Load Telogator allele data
    df = pd.read_csv(
        input_file, sep="\t", header=0,
        usecols=[0, 4, 5, 8, 10],
        names=["chr", "TL_p75", "read_TLs", "tvr_len", "read_ids"]
    )

    # 3. Feature Engineering
    df["supporting_reads"] = df["read_ids"].fillna("").str.split(",").str.len()
    df["max_read_TL"] = df["read_TLs"].apply(max_read_tl)
    df["TL_p75"] = pd.to_numeric(df["TL_p75"], errors="coerce")
    df["tvr_len"] = pd.to_numeric(df["tvr_len"], errors="coerce")

    # Expand clustered chromosome arms
    df["chr"] = df["chr"].astype(str)
    df["arm"] = df["chr"].str.split(",")
    df = df.explode("arm")
    df["arm"] = df["arm"].str.strip()
    df["cluster_arm_count"] = df["chr"].str.count(",") + 1

    # ==========================================
    # 4. RANKING LOGIC (Low-Stringency DP5 Variant)
    # ==========================================
    # Identify alleles meeting the low-stringency threshold of at least 5 reads
    df["has_min_support"] = df["supporting_reads"] >= 5

    # Multi-tier sorting to find the representative allele:
    # Tier 1: Prefer alleles with >= 5 reads (True comes before False)
    # Tier 2: Prefer single chromosomes over clustered proxy chromosomes (1 before 2+)
    # Tier 3: Break ties by taking the longest predicted telomere length (p75)
    df = df.sort_values(
        ["arm", "has_min_support", "cluster_arm_count", "TL_p75"],
        ascending=[True, False, True, False]
    )

    # Assign rank based on this sorting order
    df["allele_rank"] = df.groupby("arm").cumcount() + 1
    # ==========================================

    # 5. Apply allele mode filtering (Collapsed vs Haplotype)
    if allele_mode == "collapsed":
        df = df[df["allele_rank"] == 1]
    elif allele_mode == "haplotype":
        df = df[df["allele_rank"] <= 2]

    # Flag whether the data came from a clean single chromosome or a cluster
    df["source_type"] = np.where(df["cluster_arm_count"] == 1, "single", "proxy_from_cluster")

    # Rename columns for final output formatting
    df = df.rename(columns={
        "arm": "chr_arm",
        "TL_p75": "telomere_length_based_on_TL_p75",
        "max_read_TL": "max_individual_read_TL",
        "source_type": "is_value_from_single_chr_or_clustered_chr"
    })

    out_cols = [
        "chr_arm", "allele_rank", "telomere_length_based_on_TL_p75",
        "max_individual_read_TL", "tvr_len", 
        "is_value_from_single_chr_or_clustered_chr", "cluster_arm_count"
    ]
    out = df[out_cols]

    if allele_mode == "collapsed":
        out = out.drop(columns="allele_rank")
    else:
        out = out.rename(columns={"allele_rank": "haplotype"})

    # 6. Merge with base chromosome structure to include N/As for missing arms
    base = pd.DataFrame({"chr_arm": arms})
    base["chr_gb_acc"] = base["chr_arm"].map(arm_to_gb)
    base["chr_rs_acc"] = base["chr_arm"].map(arm_to_rs)

    if allele_mode == "haplotype":
        base = base.loc[base.index.repeat(2)]
        base["haplotype"] = [1, 2] * (len(base) // 2)

    final = base.merge(
        out,
        on=["chr_arm", "haplotype"] if allele_mode == "haplotype" else ["chr_arm"],
        how="left"
    )

    final.to_csv(output_file, sep="\t", index=False, na_rep="n/a")

def main():
    parser = argparse.ArgumentParser(description="Low stringency (DP>=5) Telogator summary stats")
    parser.add_argument("--input", required=True, help="Raw tlens_by_allele.tsv file")
    parser.add_argument("--species_table", required=True, help="Reference mapping table")
    parser.add_argument("--species", required=True, help="Target sample or species ID")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--error_file", required=True, help="Error log output path")
    parser.add_argument("--allele_mode", choices=["collapsed", "haplotype"], default="collapsed",
                        help="Output single merged allele or top two haplotypes")
    args = parser.parse_args()

    process_telomere_data(args.input, args.species_table, args.species, args.output, args.error_file, args.allele_mode)

if __name__ == "__main__":
    main()
