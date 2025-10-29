#!/usr/bin/env python3
def get_options():
    import argparse
    description = "Correct ply allele based on amino acid mutation profiles"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--aa-changes", required=True, help="TSV file containing amino acid changes (aa_changes.tsv)")
    parser.add_argument("--out", required=True, help="Output TSV file with corrected ply allele numbers")
    return parser.parse_args()


if __name__ == "__main__":
    import pandas as pd
    import numpy as np
    import re

    options = get_options()

    # Load amino acid changes table
    aa_changes = pd.read_csv(options.aa_changes, sep="\t")

    # --- Detect amino acid position columns ---
    aa_cols = [c for c in aa_changes.columns if re.match(r"pos_\d+", c)]
    if not aa_cols:
        raise ValueError("No amino acid position columns detected (expected names like pos_1, pos_2, ...)")

    # --- Separate reference and sample rows ---
    ref_df = aa_changes[aa_changes["sample_id"].str.match(r"^ply-\d+$")].copy()
    sample_df = aa_changes[~aa_changes["sample_id"].str.match(r"^ply-\d+$")].copy()

    # --- Store reference patterns as tuples ---
    ref_patterns = {}
    print("\n📘 Reference Patterns:")
    for _, ref_row in ref_df.iterrows():
        ref_id = ref_row["sample_id"]
        ref_pattern = tuple(ref_row[aa_cols].fillna("."))
        ref_patterns[ref_id] = ref_pattern
        pattern_str = " | ".join([f"{c}:{v}" for c, v in zip(aa_cols, ref_pattern)])
        print(f"{ref_id}: {pattern_str}")

    # --- Add a column for ply_allele_aa if not present ---
    if "ply_allele_aa" not in aa_changes.columns:
        aa_changes.insert(aa_changes.columns.get_loc("sample_id") + 1, "ply_allele_aa", "")

    # --- Annotate each sample ---
    print("\n🔬 Matching Samples to Reference Patterns:\n")
    for idx, sample_row in sample_df.iterrows():
        sample_id = sample_row["sample_id"]
        sample_pattern = tuple(sample_row[aa_cols].fillna("."))

        # If all dots → no AA changes → ply-1
        if all(v == "." for v in sample_pattern):
            aa_changes.at[idx, "ply_allele_aa"] = "ply-1"
            print(f"Sample: {sample_id}")
            print("No AA changes → assigned ply-1\n")
            continue

        # Compare against stored reference patterns
        matched_ref = None
        for ref_id, ref_pattern in ref_patterns.items():
            if sample_pattern == ref_pattern:
                matched_ref = ref_id
                break

        # Print results
        print(f"Sample: {sample_id}")
        print("Sample pattern: ", " | ".join([f"{c}:{v}" for c, v in zip(aa_cols, sample_pattern)]))
        if matched_ref:
            print("Matched reference:", matched_ref)
            print("Ref pattern:     ", " | ".join([f"{c}:{v}" for c, v in zip(aa_cols, ref_patterns[matched_ref])]))
            aa_changes.at[idx, "ply_allele_aa"] = matched_ref
        else:
            print("No exact match found → left empty")
            aa_changes.at[idx, "ply_allele_aa"] = ""
        print()

    # --- Save corrected output ---
    aa_changes.to_csv(options.out, sep="\t", index=False)
    print(f"\n✅ Annotated table saved to: {options.out}")
