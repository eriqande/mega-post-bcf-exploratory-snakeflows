


rule angsd_do_asso_single:
	input: 
		beag="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/sections/{scaff_grp}-beagle-post.gz",
		fai=".test/config/angsd-names.fasta.fai",
		sampleFile=".test/config/dot_samples3.tsv"
	params:
		sg="{scaff_grp}",
		whichCov="sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12",
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