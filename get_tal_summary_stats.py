#!/usr/bin/env python3

"""
Telomere Length (TL) Summary Statistics Generator
-------------------------------------------------
This script processes raw Telogator allele-specific telomere length outputs.
It maps the lengths to specific chromosome arms (p and q) using a reference table,
filters for high-quality reads, and selects the top alleles based on the chosen 
mode ('collapsed' for a single representative allele, or 'haplotype' for two).
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
    """Parses a comma-separated string of read lengths and returns the maximum value."""
    try:
        vals = [int(x) for x in str(tl_string).split(",") if x.strip()]
        return max(vals) if vals else np.nan
    except Exception:
        return np.nan

def load_chr_arms(species_table, species):
    """
    Loads the chromosome mapping table and generates expected 'p' and 'q' arms 
    for the specified species/sample. Returns mapping dictionaries for accession IDs.
    """
    df = pd.read_csv(species_table, sep="\t")
    
    # Filter for the exact species/sample ID
    df = df[df["REF_ACC"] == species]
    if df.empty:
        sys.exit(f"Species not found in chromosome table: {species}")

    arms = []
    arm_to_gb = {}
    arm_to_rs = {}

    for _, r in df.iterrows():
        chrom = r["CHROM"]
        # Create both p and q arm entries for each chromosome
        for suffix in ["p", "q"]:
            arm = f"{chrom}{suffix}"
            arms.append(arm)
            arm_to_gb[arm] = r["GB_ACC"]
            arm_to_rs[arm] = r["RS_ACC"]

    return natural_sort(arms), arm_to_gb, arm_to_rs

def process_telomere_data(input_file, species_table, species, output_file, error_file, allele_mode="collapsed"):
    """Core function to read, filter, rank, and format telomere length data."""
    arms, arm_to_gb, arm_to_rs = load_chr_arms(species_table, species)

    # 1. Handle missing input files gracefully by generating an empty formatted table
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

    # 2. Load the specific columns needed from the Telogator TSV
    df = pd.read_csv(
        input_file, sep="\t", header=0,
        usecols=[0, 4, 5, 8, 10],
        names=["chr", "TL_p75", "read_TLs", "tvr_len", "read_ids"]
    )

    # 3. Feature Engineering: Calculate read support and maximum individual read lengths
    df["supporting_reads"] = df["read_ids"].fillna("").str.split(",").str.len()
    df["max_read_TL"] = df["read_TLs"].apply(max_read_tl)
    df["TL_p75"] = pd.to_numeric(df["TL_p75"], errors="coerce")
    
    # Expand clustered chromosomes (e.g., 'chr1,chr2' becomes two rows)
    df["chr"] = df["chr"].astype(str)
    df["arm"] = df["chr"].str.split(",")
    df = df.explode("arm")
    df["arm"] = df["arm"].str.strip()
    df["cluster_arm_count"] = df["chr"].str.count(",") + 1

    # 4. Rank alleles based on read support (highest support ranks first)
    df = df.sort_values(["arm", "supporting_reads"], ascending=[True, False])
    df["allele_rank"] = df.groupby("arm").cumcount() + 1

    # 5. Filter based on user-selected mode
    if allele_mode == "collapsed":
        df = df[df["allele_rank"] == 1]
    elif allele_mode == "haplotype":
        df = df[df["allele_rank"] <= 2]

    df["source_type"] = np.where(df["cluster_arm_count"] == 1, "single", "proxy_from_cluster")

    # 6. Merge with the base chromosome structure to ensure missing arms are represented as N/A
    base = pd.DataFrame({"chr_arm": arms})
    if allele_mode == "haplotype":
        base = base.loc[base.index.repeat(2)]
        base["haplotype"] = [1, 2] * (len(base) // 2)

    df = df.rename(columns={"arm": "chr_arm"})
    final = base.merge(
        df, 
        on=["chr_arm", "allele_rank"] if allele_mode == "haplotype" else ["chr_arm"], 
        how="left"
    )

    final.to_csv(output_file, sep="\t", index=False, na_rep="n/a")

def main():
    parser = argparse.ArgumentParser(description="Summarize Telogator allele outputs")
    parser.add_argument("--input", required=True, help="Raw tlens_by_allele.tsv file")
    parser.add_argument("--species_table", required=True, help="Reference mapping table")
    parser.add_argument("--species", required=True, help="Target sample or species ID")
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("--error_file", required=True, help="Error log output path")
    parser.add_argument("--allele_mode", choices=["collapsed", "haplotype"], default="collapsed", 
                        help="Output a single merged allele or top two haplotypes")
    args = parser.parse_args()

    process_telomere_data(args.input, args.species_table, args.species, args.output, args.error_file, args.allele_mode)

if __name__ == "__main__":
    main()
