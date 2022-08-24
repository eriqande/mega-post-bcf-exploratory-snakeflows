import os
import warnings
import pandas as pd
from snakemake.utils import validate



# get all the scaff_group_ids
scaffold_groups = pd.read_table(config["scaff_group_path"]).set_index("id", drop=False)
# ensure that column order is correct
scaff_cols = list(scaffold_groups.columns)
if scaff_cols[0] != 'id' or scaff_cols[1] != 'chrom' or scaff_cols[2] != 'start' or scaff_cols[3] != 'stop' or scaff_cols[4] != "angsd_chrom": 
	raise Exception("Column order is important in the scaffold_groups file.  The columns must be id, chrom, start, stop, in that order.")
unique_scaff_groups=list(scaffold_groups.id.unique())

# then, also get the first scaff group id
first_scaff_group_id=unique_scaff_groups[0]



def get_bcf_path(wildcards):
	return config["bcf"][wildcards.bcf_id]["path"]

def get_subsamp_path(wildcards):
	return config["bcf"][wildcards.bcf_id]["sample_subsets"][wildcards.sampsub]["path"]

def get_scaff_group_path(wildcards):
	return config["scaff_group_path"]

def get_bcftools_opts(wildcards):
	return config["bcftools_opts"][wildcards.bcfilt]
