#!/usr/bin/env python3

import pandas as pd
import sys
import re

import argparse

parser = argparse.ArgumentParser(
    description="Extract numeric and relative abundance columns from a Bracken combine results table."
)
parser.add_argument("--input", required=True, help="Bracken combined abundance txt (tab-separated)")
parser.add_argument("--level", required=True, choices=["species", "genus", "phylum"], help="Taxonomic level")
parser.add_argument("--out-raw", required=True, help="Output: CSV with raw abundance values")
parser.add_argument("--out-rel", required=True, help="Output: CSV with relative abundance (fraction) values")
args = parser.parse_args()

df = pd.read_csv(args.input, sep='\t')

# Find the part to match for column selection
raw_suffix = f"_bracken.{args.level}.report.txt_num"
rel_suffix = f"_bracken.{args.level}.report.txt_frac"

# 1. Extract columns: first 3 + all *_num
num_columns = [c for c in df.columns[:3]] + [c for c in df.columns if c.endswith(raw_suffix)]
species_abundance_num = df[num_columns]
# Rename *_num suffix
species_abundance_num.columns = [
    re.sub(re.escape(raw_suffix)+'$', "", c)
    for c in species_abundance_num.columns
]

# 2. Extract columns: first 3 + all *_frac
frac_columns = [c for c in df.columns[:3]] + [c for c in df.columns if c.endswith(rel_suffix)]
species_abundance_frac = df[frac_columns]
# Rename *_frac suffix
species_abundance_frac.columns = [
    re.sub(re.escape(rel_suffix)+'$', "", c)
    for c in species_abundance_frac.columns
]

# Write to CSV (no row numbers/index, no quoting strings)
species_abundance_num.to_csv(args.out_raw, index=False, quoting=0)
species_abundance_frac.to_csv(args.out_rel, index=False, quoting=0)