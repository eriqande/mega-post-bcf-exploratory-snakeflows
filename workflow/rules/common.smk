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


def get_do_asso_param_set(wildcards):
	return config["params"]["do_asso"][wildcards.param_set]


wildcard_constraints:
	thin_int="[0-9]+",
	thin_start="[0-9]+"



# here are the functions used to convert elements in targets
# in the config into requested file names.
# here, analysis is the key under targets, and tl is the list
# of three things: [Main, SubSamp, Params]

# this gives us the main starting path for any analysis and triplet
def main_params_path(analysis, tl):
	tlists=config["targets"][analysis]
	return "results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_spec}/{anal}/maf_{min_maf}/{param_set}".format(
		bcf_id=config["main_params"][tl[0]]["bcf"],
		bcfilt=config["main_params"][tl[0]]["filt"],
		sampsub=tl[1],
		thin_spec=config["main_params"][tl[0]]["thin_spec"],
		anal=analysis,
		min_maf=config["main_params"][tl[0]]["maf"],
		param_set=tl[2])

# this gives us the main starting path for any analysis that doesn't have a param set
def beagle_path(analysis, tl):
	tlists=config["targets"][analysis]
	return "results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_spec}/beagle-gl/beagle-gl.gz".format(
		bcf_id=config["main_params"][tl[0]]["bcf"],
		bcfilt=config["main_params"][tl[0]]["filt"],
		sampsub=tl[1],
		thin_spec=config["main_params"][tl[0]]["thin_spec"],
		anal=analysis,
		min_maf=config["main_params"][tl[0]]["maf"]
		)



# this function parses the dictionaries in config["targets"] and
# expands to all the different requested ouptut files
def expand_targets():
	targ=config["targets"]
	ret = []
	if "do_asso" in targ:
		for T in targ["do_asso"]:
			mainp = main_params_path("do_asso", T)
			# get the LRT files
			ret = ret + [mainp + "/all-scaff-groups.lrt0.gz"]
			# get the manhattan plots
			ret = ret + [mainp + "/{sampsub}-{param_set}-manhattan-plot.jpg".format(sampsub=T[1], param_set=T[2])]  
	if "beagle_gl" in targ:
		for T in targ["beagle_gl"]:
			mainp = beagle_path("beagle_gl", T)
			ret = ret + [mainp]
	if "beagle_regions" in targ:
		for T in targ["beagle_regions"]:
			path="results/beagle_regions/{mpars}/{region}/phased/{sampset}.vcf.gz.csi".format(
				mpars=T[0],
				region=T[2],
				sampset=T[1]
				)
			ret = ret + [path]
	if "pcangsd_plain" in targ:
		for T in targ["pcangsd_plain"]:
			mainp = main_params_path("pcangsd_plain", T)
			ret = ret + [mainp + "/out.cov"]
	return ret













