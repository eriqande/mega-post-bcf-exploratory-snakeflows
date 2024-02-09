


# a simple rule to have angsd produce mafs for the sample subset.  
# this will do it straight from the bcf file
rule get_mafs_from_angsd:
	input:
		bcf="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/sections/{scaff_group}.bcf",
		csi="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/main.bcf.csi",
		fai=config["fai_path"]
		info="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/info/samples.txt"
	output:
		mafs="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/mafs/sections/{scaff_group}.mafs.gz"
	conda:
		"../envs/angsd.yaml"
	log:
		"results/logs/get_mafs_from_angsd/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/{scaff_group}.log"
	benchmark:
		"results/benchmarks/get_mafs_from_angsd/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/{scaff_group}.bmk"
	shell:
		" NIND=$(awk 'NF>0 {{++n}} END {{print n}}' {input.info}); "
		" PREFIX=$(echo {output.mafs} | sed 's/\.mafs\.gz//g;'); "
		" angsd -vcf-pl {input.bcf} -fai {input.fai} -nind $NIND -domaf 3 -out $PREFIX 2> {log} "
