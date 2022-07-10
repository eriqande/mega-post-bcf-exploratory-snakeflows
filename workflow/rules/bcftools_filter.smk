


rule bcf_samps_and_filt_scatter:
	input:
		bcf=get_parent_bcf_from_id,
		samps=get_sample_subset_path,
		scaff_grp_path=get_scaff_group_file_for_bcf
	params:
		bcftools_opts=get_bcftools_opts,
		sg="{scaff_grp}"
	output:
		scaff_members="results/bcf/{bcf_id}/scaff_members/{scaff_grp}.scaff_members.tsv",
		bcf=temp("results/bcf/{bcf_id}/thin_0_0/sections/{scaff_grp}.bcf"),
		stats=temp("results/bcf/{bcf_id}/thin_0_0/sections/{scaff_grp}.bcf_stats.txt"),
		pos="results/bcf/{bcf_id}/thin_0_0/sections/{scaff_grp}.positions.tsv.gz",
	conda:
		"../envs/bcftools.yaml"
	log:
		"results/logs/bcf_samps_and_filt_scatter/{bcf_id}/sections/{scaff_grp}.log"
	benchmark:
		"results/benchmarks/bcf_samps_and_filt_scatter/{bcf_id}/sections/{scaff_grp}.bmk"
	shell:
		" (awk -v sg='{params.sg}' -f workflow/scripts/get_scaff_members.awk {input.scaff_grp_path} > {output.scaff_members} && "
		"    bcftools view -Ou -R {output.scaff_members} -S {input.samps} {input.bcf} | "
		"    bcftools +fill-tags -Ou -- -t all | "
		"    bcftools view -Ob {params.bcftools_opts} > {output.bcf}  && "
		"    bcftools stats --af-tag MAF {output.bcf} > {output.stats} && "
		"    bcftools query -f '%CHROM\\t%POS\\n' {output.bcf} | gzip -c > {output.pos} "
		" ) 2> {log} "



rule bcf_samps_and_filt_gather:
	input:
		bcfs=lambda wc: expand("results/bcf/{{bcf_id}}/thin_0_0/sections/{sg}.bcf", sg=all_scaff_group_ids(wc)),
		poses=lambda wc: expand("results/bcf/{{bcf_id}}/thin_0_0/sections/{sg}.positions.tsv.gz", sg=all_scaff_group_ids(wc)),
		statses=lambda wc: expand("results/bcf/{{bcf_id}}/thin_0_0/sections/{sg}.bcf_stats.txt", sg=all_scaff_group_ids(wc)),
		samps_path=get_sample_subset_path,
		scaff_grp_path=get_scaff_group_file_for_bcf,
		bcf_tsv=config["bcf_tsv"]
	output:
		bcf="results/bcf/{bcf_id}/thin_0_0/main.bcf",
		csi="results/bcf/{bcf_id}/thin_0_0/main.bcf.csi",
		info="results/bcf/{bcf_id}/thin_0_0/info/info.txt",
		pos="results/bcf/{bcf_id}/thin_0_0/info/positions.tsv.gz",
		stats="results/bcf/{bcf_id}/thin_0_0/info/bcf_stats.txt",
		samp_sub_file="results/bcf/{bcf_id}/thin_0_0/info/sample_subset_file.txt",
		scaff_grp_file="results/bcf/{bcf_id}/thin_0_0/info/scaffold_groups.tsv",
		samps="results/bcf/{bcf_id}/thin_0_0/info/samples.txt"
	conda:
		"../envs/bcftools.yaml"
	log:
		"results/logs/bcf_samps_and_filt_gather/{bcf_id}/main.log"
	benchmark:
		"results/benchmarks/bcf_samps_and_filt_gather/{bcf_id}/main.bmk"
	shell:
		"( bcftools concat --naive {input.bcfs} > {output.bcf} && "
		"  bcftools index {output.bcf} && "
		"  cat {input.poses} > {output.pos} && "
		"  plot-vcfstats -m {input.statses} > {output.stats}  && "
		"  cp {input.samps_path} {output.samp_sub_file}  && "
		"  cp {input.scaff_grp_path} {output.scaff_grp_file}  && "
		"  bcftools query -l {output.bcf} > {output.samps}  && "
		"  awk -v id='{wildcards.bcf_id}' -f workflow/scripts/tsv_row2yaml.awk {input.bcf_tsv} > {output.info} "
		") 2> {log} "

