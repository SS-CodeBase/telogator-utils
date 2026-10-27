#!/usr/bin/env python3

"""
Autosomal Telomere Length Aggregator (Multi-Condition)
------------------------------------------------------
This script compiles the per-sample chromosome summaries generated 
in Step 1 into a single master matrix. It filters out sex chromosomes 
(X and Y) to prevent sex-bias in the cohort-level data, and computes 
standard statistical metrics (Min, Max, Median, Mean, StDev) across 
five different combinations of Telomere and TVR lengths.
"""

import pandas as pd
import numpy as np
import argparse
import glob
import os
import re

def is_autosome(chr_arm):
    """
    Checks if a chromosome arm is an autosome.
    Uses regex to filter out sex chromosomes like chrX, chrY, chrXp, chrYq.
    Returns True if autosomal, False otherwise.
    """
    return not re.match(r"chr[XY]\d*[pq]$", chr_arm, re.IGNORECASE)

def compute_stats(values):
    """
    Computes summary statistics for a given pandas Series of numeric values.
    Returns a list: [Min, Max, Median, Average, Standard Deviation].
    Requires at least 2 values to compute standard deviation.
    """
    if len(values) == 0:
        return ["n/a"] * 5

    return [
        round(values.min(), 2),
        round(values.max(), 2),
        round(values.median(), 2),
        round(values.mean(), 2),
        round(values.std(ddof=1), 2) if len(values) > 1 else "n/a"
    ]

def main():
    parser = argparse.ArgumentParser(description="Aggregate autosomal telomere stats across samples")
    parser.add_argument("--input_dir", required=True, help="Directory containing Step 1 TSV outputs")
    parser.add_argument("--output", required=True, help="Path for the final aggregated TSV matrix")
    args = parser.parse_args()

    # Find all summary TSVs matching the expected output format of Step 1
    summary_files = glob.glob(os.path.join(args.input_dir, "*.all.chr.tal.summary.*.tsv"))

    results = []
    no_tal_samples = []

    # =========================================================================
    # Process Each Sample File
    # =========================================================================
    for file in sorted(summary_files):
        # Extract the sample name from the filename
        sample = re.sub(r"\.all\.chr\.tal\.summary\..*\.tsv$", "", os.path.basename(file))

        try:
            df = pd.read_csv(file, sep="\t")

            if df.empty:
                no_tal_samples.append(sample)
                continue

            # Isolate autosomes only
            autosomes = df[df["chr_arm"].apply(is_autosome)].copy()
            # Uncomment the next line if you need to strictly avoid haplotype double-counting:
            # autosomes = autosomes.drop_duplicates(subset=["chr_arm"]) 

            if autosomes.empty:
                no_tal_samples.append(sample)
                continue

            # Ensure core metrics are treated as numeric for mathematical operations
            numeric_cols = ["telomere_length_based_on_TL_p75", "tvr_len", "max_individual_read_TL"]
            for col in numeric_cols:
                if col in autosomes.columns:
                    autosomes[col] = pd.to_numeric(autosomes[col], errors="coerce")

            # Filter out N/A and invalid negative predictions
            autosomes = autosomes[autosomes["telomere_length_based_on_TL_p75"].notna()]
            autosomes = autosomes[autosomes["telomere_length_based_on_TL_p75"] >= 0]

            if autosomes.empty:
                no_tal_samples.append(sample)
                continue

            # Compute general dataset health metrics for this sample
            total_autosomal_arms = len(df[df["chr_arm"].apply(is_autosome)])
            tal_count = len(autosomes)
            percentage = (tal_count / total_autosomal_arms) * 100 if total_autosomal_arms > 0 else 0

            # =================================================================
            # Define the 5 Analytical Conditions
            # =================================================================
            cond1 = autosomes["telomere_length_based_on_TL_p75"] + autosomes["tvr_len"].fillna(0)
            cond2 = autosomes["telomere_length_based_on_TL_p75"]
            cond3 = autosomes["max_individual_read_TL"]
            cond4 = autosomes["max_individual_read_TL"] + autosomes["tvr_len"].fillna(0)
            cond5 = autosomes["tvr_len"]

            # Initialize the row data with base metadata
            row = [sample, tal_count, total_autosomal_arms, round(percentage, 2)]

            # Compute and append statistics for all 5 conditions sequentially
            for cond in [cond1, cond2, cond3, cond4, cond5]:
                cond = cond.dropna()
                row.extend(compute_stats(cond))

            results.append(row)

        except Exception as e:
            print(f"Skipping {sample} due to error: {e}")
            no_tal_samples.append(sample)

    # =========================================================================
    # Format and Save Output
    # =========================================================================
    # Base columns
    columns = [
        "Sample",
        "Number_of_autosome_chr_arms_TAL_prediction",
        "Total_no_of_autosome_arms",
        "Percentage_of_chr_ends_with_prediction",
    ]

    # Generate column headers dynamically for the 5 conditions
    for i in range(1, 6):
        columns.extend([
            f"Cond{i}_Minimum",
            f"Cond{i}_Maximum",
            f"Cond{i}_Median",
            f"Cond{i}_Average",
            f"Cond{i}_Standard_Deviation"
        ])

    out_df = pd.DataFrame(results, columns=columns)
    out_df.sort_values("Sample", inplace=True)
    
    # Write the master matrix
    out_df.to_csv(args.output, sep="\t", index=False)

    # Write the missing samples log file to the same directory as the output
    log_file = os.path.join(os.path.dirname(args.output), "no_tal_values.samples.log")
    with open(log_file, "w") as f:
        for sample in sorted(no_tal_samples):
            f.write(sample + "\n")

    print(f"Aggregation complete → {args.output}")
    print(f"Samples without autosomal TAL values logged → {log_file}")

if __name__ == "__main__":
    main()
