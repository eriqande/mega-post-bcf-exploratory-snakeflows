import os
import warnings
import pandas as pd
from snakemake.utils import validate


bcf_table = pd.read_table(config["bcf_tsv"]).set_index("id", drop=False)
validate(bcf_table, schema="../schemas/bcf_tsv.schema.yaml")


# get a list of just the unique values of the scaffold_group and of the chromosomes
#unique_scaff_groups = list(scaffold_groups.id.unique())


bcftools_opts_table = pd.read_table(config["bcftools_opts"]).set_index("id", drop=False)
validate(bcftools_opts_table, schema="../schemas/bcftools_opts.schema.yaml")


# get bcf path from the bcf_file wildcard
def bcf_from_ID(wildcards):
	return bcf_table.loc[ wildcards.bcf_file, "bcf_path" ]

def bcf_csi_from_ID(wildcards):
	return "{path}.csi".format(path = bcf_table.loc[ wildcards.bcf_file, "bcf_path" ])

# get the path to the sample subset file from the sample_subset
def sample_subset_file(wildcards):
	return "{dir}/{ss}.txt".format(dir = config['samp_subs_dir'], ss = wildcards.sample_subset)

# get the bcftools options from the bcftools_opts wildcard
def get_bcftools_opts(wildcards):
	return bcftools_opts_table.loc[ wildcards.bcftools_opts, "opts" ]


# get scaffold groups file path for a specific BCF file
def get_scaff_group_file_for_bcf(wildcards):
	return bcf_table.loc[ wildcards.bcf_file, "scaff_group_path" ]

# return all the scaffold groups given a bcf_file id
def return_scaffold_groups_table(wildcards):
	sg_file=get_scaff_group_file_for_bcf(wildcards)
	scaffold_groups = pd.read_table(sg_file).set_index("id", drop=False)
	#validate(scaffold_groups, schema="../schemas/scaffold_groups.schema.yaml")
	# ensure that column 1 of the scaffold group file is "id" and
	# column 2 is "chrom".  This is essential because I used those
	# column positions in some awk to pull things out.
	scaff_cols = list(scaffold_groups.columns)
	if scaff_cols[0] != 'id' or scaff_cols[1] != 'chrom' or scaff_cols[2] != 'start' or scaff_cols[3] != 'stop': 
		raise Exception("Column order is important in the scaffold_groups file.  The columns must be id, chrom, start, stop, in that order.")
	return(scaffold_groups)

# return just the first possible scaff group, given the bcf_file id
def first_scaff_group_id(wildcards):
	stable=return_scaffold_groups_table(wildcards)
	return stable.iloc[0].id

# all the unique scaff group ids, in order
def all_scaff_group_ids(wildcards):
	stable=return_scaffold_groups_table(wildcards)
	return list(stable.id.unique())
