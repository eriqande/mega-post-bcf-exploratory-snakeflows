# Rules related to switching between formats and things




rule beagle3_glikes_from_PL_scattered:
	input:
		bcf=bcf_from_ID,
		csi=bcf_csi_from_ID,
		regions="results/bcf_{bcf_file}/bcf_region_files/{scaff_grp}.tsv"
	params:
		sfile=sample_subset_file,
		bcfopts=get_bcftools_opts
	output:
		top_row="results/bcf_{bcf_file}/thin_{thin_int}_{thin_start}/samp-sub_{sample_subset}/bcft-opts_{bcftools_opts}/beagle-gl/sections/{scaff_grp}.toprow.gz",
		body="results/bcf_{bcf_file}/thin_{thin_int}_{thin_start}/samp-sub_{sample_subset}/bcft-opts_{bcftools_opts}/beagle-gl/sections/{scaff_grp}.body.gz"
	log:
		"results/logs/beagle3_glikes_from_PL_scattered/bcf_{bcf_file}/thin_{thin_int}_{thin_start}/samp-sub_{sample_subset}/bcft-opts_{bcftools_opts}/sections/{scaff_grp}.log"
	conda:
		"../envs/bcftools.yaml"
	shell:
		" (bcftools query -S {params.sfile} -l {input.bcf} | " 
		" awk -f workflow/scripts/beagle3header.awk  | gzip -c > {output.top_row}) 2> {log}; "
		"        "
		" (bcftools view -Ou {params.bcfopts} -S {params.sfile} -R {input.regions} {input.bcf} |  "
		" bcftools query -f '%CHROM:%POS\\t%REF\\t%ALT[\\t%PL]\\n' | "
		" awk -f workflow/scripts/pl2gl.awk | gzip -c  >  {output.body}) 2>> {log} "



rule gather_beagle3_glikes_from_PL:
	input: 
		header=lambda wc: expand("results/bcf_{{bcf_file}}/thin_{{thin_int}}_{{thin_start}}/samp-sub_{{sample_subset}}/bcft-opts_{{bcftools_opts}}/beagle-gl/sections/{sg}.toprow.gz", sg=first_scaff_group_id(wc)),
		scaff_gzs = lambda wc: expand("results/bcf_{{bcf_file}}/thin_{{thin_int}}_{{thin_start}}/samp-sub_{{sample_subset}}/bcft-opts_{{bcftools_opts}}/beagle-gl/sections/{sg}.body.gz", sg=all_scaff_group_ids(wc))
	output:
		"results/bcf_{bcf_file}/thin_{thin_int}_{thin_start}/samp-sub_{sample_subset}/bcft-opts_{bcftools_opts}/beagle-gl/all-scaffs.gz"
	shell:
		"cat {input.header} {input.scaff_gzs} > {output}"