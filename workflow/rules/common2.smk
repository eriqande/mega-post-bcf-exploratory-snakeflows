import os
import warnings
import pandas as pd
from snakemake.utils import validate




def get_bcf_path(wildcards):
	return config["bcf"][wildcards.bcf_id]["path"]

def get_subsamp_path(wildcards):
	return config["bcf"][wildcards.bcf_id]["sample_subsets"][wildcards.sampsub]["path"]

def get_scaff_group_path(wildcards):
	return config["bcf"][wildcards.bcf_id]["scaff_group_path"]

def get_bcftools_opts(wildcards):
	return config["bcftools_opts"][wildcards.bcfilt]

def all_scaff_group_ids(wildcards):
	sgpath=config["bcf"][wildcards.bcf_id]["scaff_group_path"]
	scaffold_groups = pd.read_table(sgpath).set_index("id", drop=False)
	# ensure that column order is correct
	scaff_cols = list(scaffold_groups.columns)
	if scaff_cols[0] != 'id' or scaff_cols[1] != 'chrom' or scaff_cols[2] != 'start' or scaff_cols[3] != 'stop' or scaff_cols[4] != "angsd_chrom": 
		raise Exception("Column order is important in the scaffold_groups file.  The columns must be id, chrom, start, stop, in that order.")
	return list(scaffold_groups.id.unique())

def first_scaff_group_id(wildcards):
	all_them=all_scaff_group_ids(wildcards)
	return all_them[0]


# this is a hack that remains so I can make output files as needed.
# I've gotta get a better way around this, or hardwire the scaff groups
# file
stable = pd.read_table(".test/bcf/test.bcf.scaff_groups.tsv").set_index("id", drop=False)
unique_scaff_groups = list(stable.id.unique())