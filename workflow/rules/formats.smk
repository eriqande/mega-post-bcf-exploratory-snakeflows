# Rules related to switching between formats and things




# makes beagle GL file from the PL field in BCF file
rule bcf2beagle_gl_scatter:
	input:
		bcf="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/main.bcf",
		csi="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/main.bcf.csi",
		regions="results/bcf/{bcf_id}/scaff_members/{scaff_grp}.scaff_members.tsv",
		sfile="results/bcf/{bcf_id}/thin_{thin_int}_{thin_start}/info/samples.txt"
	output:
		body=temp("results/beagle-gl/{bcf_id}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.body.gz"),
		top_row=temp("results/beagle-gl/{bcf_id}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.toprow.gz"),
	log:
		"results/logs/bcf2beagle_gl_scatter/bcf_{bcf_id}/thin_{thin_int}_{thin_start}/{scaff_grp}.log"
	benchmark:
		"results/benchmarks/bcf2beagle_gl_scatter/bcf_{bcf_id}/thin_{thin_int}_{thin_start}/{scaff_grp}.bmk"
	conda:
		"../envs/bcftools.yaml"
	shell:
		" (           " 
		"    awk -f workflow/scripts/beagle3header.awk {input.sfile} | gzip -c > {output.top_row}  && "
		"    bcftools view -Ou -R {input.regions} {input.bcf} |  "
		"    bcftools query -f '%CHROM:%POS\\t%REF\\t%ALT[\\t%PL]\\n' | "
		"    awk -f workflow/scripts/pl2gl.awk | gzip -c  >  {output.body} " 
		" ) 2> {log}  "



rule bcf2beagle_gl_gather:
	input: 
		header=lambda wc: expand("results/beagle-gl/{{bcf_id}}/thin_{{thin_int}}_{{thin_start}}/sections/{sg}.toprow.gz", sg=first_scaff_group_id(wc)),
		scaff_gzs = lambda wc: expand("results/beagle-gl/{{bcf_id}}/thin_{{thin_int}}_{{thin_start}}/sections/{sg}.body.gz", sg=all_scaff_group_ids(wc))
	output:
		"results/beagle-gl/{bcf_id}/thin_{thin_int}_{thin_start}/beagle-gl.gz"
	log:
		"results/logs/bcf2beagle_gl_gather/bcf_{bcf_id}/thin_{thin_int}_{thin_start}/main.log"
	benchmark:
		"results/logs/bcf2beagle_gl_gather/bcf_{bcf_id}/thin_{thin_int}_{thin_start}/main.bmk"
	shell:
		"cat {input.header} {input.scaff_gzs} > {output} 2> {log} "




# this one pastes the gzipped top row and gzipped body together
# to make a whole beagle gl file.  It is putting a header on it, so
# we call it capitating.
rule capitate_beagle_gl_sections:
	input: 
		body="results/beagle-gl/{bcf_id}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.body.gz",
		top_row="results/beagle-gl/{bcf_id}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.toprow.gz",
	output:
		beag="results/beagle-gl/{bcf_id}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.beagle-gl.gz",
	log:
		"results/logs/capitate_beagle_gl_sections/bcf_{bcf_id}/thin_{thin_int}_{thin_start}/{scaff_grp}.log"
	benchmark:
		"results/benchmarks/capitate_beagle_gl_sections/bcf_{bcf_id}/thin_{thin_int}_{thin_start}/{scaff_grp}.bmk"
	shell:
		"cat {input.top_row} {input.body} > {output.beag} 2> {log} "

