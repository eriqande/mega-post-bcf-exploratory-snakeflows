

rule beagle_regions_bcftools_filter_and_subsamp:
	input:
		bcf = lambda wildcards: config["bcf"][config["main_params"][wildcards.br_main]["bcf"]]["path"],
		subsamp = lambda wildcards: config["bcf"][config["main_params"][wildcards.br_main]["bcf"]]["sample_subsets"][wildcards.br_subsamp]["path"],
	params:
		region = lambda wildcards: config["params"]["beagle_regions"]["regions"][wildcards.br_reg],
		filt = lambda wildcards: config["bcftools_opts"][config["main_params"][wildcards.br_main]["filt"]]
	output:
		vcf="results/beagle_regions/{br_main}/{br_reg}/vcf/{br_subsamp}.vcf.gz",
	conda:
		"../envs/bcftools.yaml"
	log:
		"results/logs/beagle_regions_bcftools_filter_and_subsamp/{br_main}/{br_reg}/{br_subsamp}.log"
	benchmark:
		"results/benchmarks/beagle_regions_bcftools_filter_and_subsamp/{br_main}/{br_reg}/{br_subsamp}.log"
	shell:
		"(bcftools view -S {input.subsamp} -Ou {input.bcf} {params.region} | "
		"  bcftools view -Oz {params.filt} > {output.vcf}) 2> {log} "


rule impute_em_beagle4:
	input:
		vcf="results/beagle_regions/{br_main}/{br_reg}/vcf/{br_subsamp}.vcf.gz",
	params:
		region=lambda wildcards: config["params"]["beagle_regions"]["regions"][wildcards.br_reg] 
	output:
		imp="results/beagle_regions/{br_main}/{br_reg}/imputed/{br_subsamp}.vcf.gz",
	conda:
		"../envs/beagle41.yaml"
	log:
		"results/logs/impute_em_beagle4/{br_main}/{br_reg}/{br_subsamp}.log"
	benchmark:
		"results/benchmarks/impute_em_beagle4/{br_main}/{br_reg}/{br_subsamp}.log"
	shell:
		" OUTPRE=$(echo {output.imp} | sed 's/.vcf.gz$//g;');     "
		" beagle gl={input.vcf} out=$OUTPRE > {log} 2>&1 "




rule phase_em_beagle4:
	input:
		imp="results/beagle_regions/{br_main}/{br_reg}/imputed/{br_subsamp}.vcf.gz",
	params:
		region= lambda wildcards: config["params"]["beagle_regions"]["regions"][wildcards.br_reg] 
	output:
		vcf="results/beagle_regions/{br_main}/{br_reg}/phased/{br_subsamp}.vcf.gz",
	conda:
		"../envs/beagle41.yaml"
	log:
		"results/logs/phase_em_beagle4/{br_main}/{br_reg}/{br_subsamp}.log"
	benchmark:
		"results/benchmarks/phase_em_beagle4/{br_main}/{br_reg}/{br_subsamp}.log"
	shell:
		" OUTPRE=$(echo {output} | sed 's/.vcf.gz$//g;');     "
		" beagle gt={input.imp} out=$OUTPRE > {log} 2>&1 "



rule index_phased_vcfs:
	input:
		vcf="results/beagle_regions/{br_main}/{br_reg}/phased/{br_subsamp}.vcf.gz"
	output:
		idx="results/beagle_regions/{br_main}/{br_reg}/phased/{br_subsamp}.vcf.gz.csi"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools index {input} "
		""