#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_options():
    parser = argparse.ArgumentParser(description="Update FASTA headers based on mapping file (4th column)")
    parser.add_argument("--input-fasta", required=True, help="Input FASTA file")
    parser.add_argument("--mapping-file", required=True, help="Mapping file (tab-delimited, 4th column used)")
    parser.add_argument("--output-fasta", required=True, help="Output FASTA file with updated headers")
    return parser.parse_args()

def load_mapping(mapping_file):
    """
    Load the mapping file. Returns a dict mapping original prefix to ply number.
    """
    mapping = {}
    with open(mapping_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                sample_col = parts[0]  # e.g., '680_ply-2'
                ply_col = parts[1]     # e.g., 'ply-2'

                # Extract numeric sample ID
                if "_ply-" in sample_col:
                    sample_id = sample_col.split("_ply-")[0]
                else:
                    sample_id = sample_col

                # Extract ply number from ply_col
                if "ply-" in ply_col:
                    ply_number = ply_col.split("ply-")[1]
                else:
                    ply_number = ply_col

                mapping[sample_id] = f"ply-{ply_number}"
    return mapping

def update_fasta_headers(input_fasta, mapping, output_fasta):
    updated_records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Extract sample ID from FASTA header (assuming header is like '680_ply-2')
        if "_" in record.id and "ply-" in record.id:
            sample_id = record.id.split("_ply-")[0]
            if sample_id in mapping:
                record.id = f"{sample_id}_{mapping[sample_id]}"
                record.description = ""
        updated_records.append(record)
    SeqIO.write(updated_records, output_fasta, "fasta")

if __name__ == "__main__":
    options = get_options()
    mapping = load_mapping(options.mapping_file)
    update_fasta_headers(options.input_fasta, mapping, options.output_fasta)
