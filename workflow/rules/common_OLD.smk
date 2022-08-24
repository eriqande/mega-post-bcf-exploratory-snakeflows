import os
import warnings
import pandas as pd
from snakemake.utils import validate



bcf_table = pd.read_table(config["bcf_tsv"]).set_index("id", drop=False)
validate(bcf_table, schema="../schemas/bcf_tsv.schema.yaml")

# get the bcftools options from the bcftools_opts column of the bcf_table
def get_bcftools_opts(wildcards):
	return bcf_table.loc[ wildcards.bcf_id, "bcftools_opts" ]

# get the parent bcf path from the bcf_id wildcard
def get_parent_bcf_from_id(wildcards):
	return bcf_table.loc[ wildcards.bcf_id, "parent_bcf_path" ]

# get the path to the sample subsets file
def get_sample_subset_path(wildcards):
	return bcf_table.loc[ wildcards.bcf_id, "sample_subset_path" ]

# get a list of just the unique values of the scaffold_group and of the chromosomes
#unique_scaff_groups = list(scaffold_groups.id.unique())



def bcf_csi_from_ID(wildcards):
	return "{path}.csi".format(path = bcf_table.loc[ wildcards.bcf_file, "bcf_path" ])

# get the path to the sample subset file from the sample_subset
def sample_subset_file(wildcards):
	return "{dir}/{ss}.txt".format(dir = config['samp_subs_dir'], ss = wildcards.sample_subset)


## Now, here we have a bit of relatively complicated stuff so
## that we can get the members of the scaffold groups on the fly
## by reading them out of a file that is determined by the wildcard.
## Though I may not need this, cuz I can use awk.

# get scaffold groups file path for a specific BCF file
def get_scaff_group_file_for_bcf(wildcards):
	return bcf_table.loc[ wildcards.bcf_id, "scaff_group_path" ]

# return all the scaffold groups given a bcf_file id
def return_scaffold_groups_table(wildcards):
	sg_file=get_scaff_group_file_for_bcf(wildcards)
	scaffold_groups = pd.read_table(sg_file).set_index("id", drop=False)
	# ensure that column order is correct
	scaff_cols = list(scaffold_groups.columns)
	if scaff_cols[0] != 'id' or scaff_cols[1] != 'chrom' or scaff_cols[2] != 'start' or scaff_cols[3] != 'stop' or scaff_cols[4] != "angsd_chrom": 
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




## Wilcard constraints.
# thin_int and thin_start must be integers greater than 0
wildcard_constraints:
	thin_int="[0-9]+",
	thin_start="[0-9]+"




## This is a total hack that I am doing to get some results for the
## Yukon Chinook before revamping all of this so that the BCF and all the associated
## files are set in the config
stable = pd.read_table(".test/bcf/test.bcf.scaff_groups.tsv").set_index("id", drop=False)
unique_scaff_groups = list(stable.id.unique())

