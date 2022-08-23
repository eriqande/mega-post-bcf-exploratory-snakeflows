


rule bcf_samps_and_filt_scatter:
	input:
		bcf=get_bcf_path,
		samps=get_subsamp_path,
		scaff_grp_path=get_scaff_group_path
	params:
		bcftools_opts=get_bcftools_opts,
		sg="{scaff_grp}"
	output:
		scaff_members="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/scaff_members/{scaff_grp}.scaff_members.tsv",
		bcf="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_0_0/sections/{scaff_grp}.bcf",
		stats="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_0_0/sections/{scaff_grp}.bcf_stats.txt",
		pos="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_0_0/sections/{scaff_grp}.positions.tsv.gz",
	conda:
		"../envs/bcftools.yaml"
	log:
		"results/logs/bcf_samps_and_filt_scatter/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/sections/{scaff_grp}.log"
	benchmark:
		"results/benchmarks/bcf_samps_and_filt_scatter/{bcf_id}/filt_{bcfilt}/{sampsub}/sections/{scaff_grp}.bmk"
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
		bcfs=lambda wc: expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_0_0/sections/{sg}.bcf", sg=all_scaff_group_ids(wc)),
		poses=lambda wc: expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_0_0/sections/{sg}.positions.tsv.gz", sg=all_scaff_group_ids(wc)),
		statses=lambda wc: expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_0_0/sections/{sg}.bcf_stats.txt", sg=all_scaff_group_ids(wc)),
		samps_path=get_subsamp_path,
		scaff_grp_path=get_scaff_group_path,
	output:
		bcf="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_0_0/main.bcf",
		csi="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_0_0/main.bcf.csi",
		pos="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_0_0/info/positions.tsv.gz",
		stats="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_0_0/info/bcf_stats.txt",
		samps="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_0_0/info/samples.txt"
	conda:
		"../envs/bcftools.yaml"
	log:
		"results/logs/bcf_samps_and_filt_gather/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/main.log"
	benchmark:
		"results/benchmarks/bcf_samps_and_filt_gather/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/main.bmk"
	shell:
		"( bcftools concat --naive {input.bcfs} > {output.bcf} && "
		"  bcftools index {output.bcf} && "
		"  cat {input.poses} > {output.pos} && "
		"  plot-vcfstats -m {input.statses} > {output.stats}  && "
		"  bcftools query -l {output.bcf} > {output.samps}  "
		") 2> {log} "

