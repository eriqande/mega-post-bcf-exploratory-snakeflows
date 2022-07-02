



rule bcftools_thin_scatter:
	input:
		bcf="results/bcf/{bcf_id}/thin_0_0/main.bcf",
		pos="results/bcf/{bcf_id}/thin_0_0/sections/{scaff_grp}.positions.tsv.gz",
		scaff_grp_path="results/bcf/{bcf_id}/thin_0_0/info/scaffold_groups.tsv"
	params:
		thin_int="{thin_int}",
		thin_start="{thin_start}",
		sg="{scaff_grp}"
	output:
		bcf=temp("results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.bcf"),
		stats=temp("results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.bcf_stats.txt"),
		pos=temp("results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.positions.tsv"),
	conda:
		"../envs/bcftools.yaml"
	log:
		"results/logs/bcftools_thin_scatter/{bcf_id}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.log"
	benchmark:
		"results/benchmarks/bcftools_thin_scatter/{bcf_id}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.bmk"
	shell:
		" (          "
		"    gunzip -c {input.pos} | awk -v thin_int={params.thin_int} -v thin_start={params.thin_start} -f workflow/scripts/thin_positions.awk > {output.pos} && "
		"    bcftools view -Ob -R {output.pos} {input.bcf}  > {output.bcf}  && "
		"    bcftools stats --af-tag MAF {output.bcf} > {output.stats} "
		" ) 2> {log} "




rule bcftools_thin_gather:
	input:
		bcfs=lambda wc: expand("results/bcf/{{bcf_id}}/thin_{{thin_int}}_{{thin_start}}/sections/{sg}.bcf", sg=all_scaff_group_ids(wc)),
		poses=lambda wc: expand("results/bcf/{{bcf_id}}/thin_{{thin_int}}_{{thin_start}}/sections/{sg}.positions.tsv", sg=all_scaff_group_ids(wc)),
		statses=lambda wc: expand("results/bcf/{{bcf_id}}/thin_{{thin_int}}_{{thin_start}}/sections/{sg}.bcf_stats.txt", sg=all_scaff_group_ids(wc)),
		samps_path="results/bcf/{bcf_id}/thin_0_0/info/samples.txt",
		scaff_grp_path="results/bcf/{bcf_id}/thin_0_0/info/scaffold_groups.tsv",
		info_file="results/bcf/{bcf_id}/thin_0_0/info/info.txt"
	output:
		bcf="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/main.bcf",
		csi="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/main.bcf.csi",
		pos="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/info/positions.tsv.gz",
		stats="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/info/bcf_stats.txt",
		samp_sub_file="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/info/sample_subset_file.txt",
		scaff_grp_file="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/info/scaffold_groups.tsv",
		samps="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/info/samples.txt",
		info_file="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/info/info.txt"
	conda:
		"../envs/bcftools.yaml"
	log:
		"results/logs/bcftools_thin_gather/{bcf_id}/thin_{thin_int}_{thin_start}/main.log"
	benchmark:
		"results/benchmarks/bcftools_thin_gather/{bcf_id}/thin_{thin_int}_{thin_start}/main.bmk"
	shell:
		"( bcftools concat --naive {input.bcfs} > {output.bcf}     && "
		"  bcftools index {output.bcf}                             && "
		"  cat {input.poses} | gzip -c  > {output.pos}             && "
		"  plot-vcfstats -m {input.statses} > {output.stats}       && "
		"  cp {input.samps_path} {output.samp_sub_file}            && "
		"  cp {input.scaff_grp_path} {output.scaff_grp_file}       && "
		"  cp {input.info_file} {output.info_file}                 && "
		"  bcftools query -l {output.bcf} > {output.samps}            "
		") 2> {log} "

