


rule angsd_do_asso_single:
	input: 
		beag="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/sections/{scaff_grp}-beagle-post.gz",
		fai=".test/config/angsd-names.fasta.fai",
		sampleFile="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/dot_samples.tsv"
	params:
		sg="{scaff_grp}",
		whichCov="cohort,PC1,PC2,PC3,PC4",
		whichPhe=" age ",
		doMaf=" -doMaf 4 ",
		what=" -doAsso 4 "
	log:
		"results/logs/do_asso/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/sections/{scaff_grp}.log"
	conda:
		"../envs/angsd.yaml"
	output:
		 arg="results/do_asso/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/sections/{scaff_grp}.arg",
		 lrt="results/do_asso/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/sections/{scaff_grp}.lrt0.gz"
	shell:
		" angsd {params.doMaf} -beagle {input.beag} -fai {input.fai} "
		"  -sampleFile {input.sampleFile} -whichPhe {params.whichPhe} "
		"  -whichCov {params.whichCov}  "
		"  -out $(dirname {output.arg})/{wildcards.scaff_grp} {params.what} > {log} 2>&1 "




rule angsd_do_asso_scatter:
	input: 
		beag="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/sections/{scaff_grp}-beagle-post.gz",
		fai="{p}-ANGSD".format(p=config["fai_path"]),
		sampleFile="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/dotsample_PCs_12.tsv"
	params:
		dicto=get_do_asso_param_set
	log:
		"results/logs/do_asso/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/maf_{min_maf}/{param_set}/sections/{scaff_grp}.log"
	conda:
		"../envs/angsd.yaml"
	output:
		 arg="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/do_asso/maf_{min_maf}/{param_set}/sections/{scaff_grp}.arg",
		 lrt="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/do_asso/maf_{min_maf}/{param_set}/sections/{scaff_grp}.lrt0.gz"
	shell:
		" angsd {params.dicto[angsd_opts]}  -beagle {input.beag} -fai {input.fai} "
		"  -sampleFile {input.sampleFile} -whichPhe {params.dicto[whichPhe]} "
		"  -whichCov {params.dicto[whichCov]}  "
		"  -out $(dirname {output.arg})/{wildcards.scaff_grp}  > {log} 2>&1 "


rule angsd_do_asso_gather:
	input:
		lrts=expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_{{thin_int}}_{{thin_start}}/do_asso/maf_{{min_maf}}/{{param_set}}/sections/{scaff_grp}.lrt0.gz", scaff_grp=unique_scaff_groups)
	output:
		lrt="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/do_asso/maf_{min_maf}/{param_set}/all-scaff-groups.lrt0.gz"
	threads: 4
	resources:
		mem_mb=19200
	shell:
		" set +o pipefail; (gunzip -c {input.lrts[0]} | head -n 1;  "
		"  for i in {input.lrts}; do gunzip -c $i | awk 'NR>1'; done "
		" ) | gzip -c > {output.lrt} "
