

# just make files that can be used with -R in bcftools 
rule write_bcf_regions:
	input:
		tsv=get_scaff_group_file_for_bcf
	params:
		sg="{scaff_grp}"
	output:
		tsv="results/bcf_{bcf_file}/bcf_region_files/{scaff_grp}.tsv"
	log:
		"results/logs/write_bcf_regions/bcf_{bcf_file}/{scaff_grp}.log"
	shell:
		" awk -F\"\\t\" 'BEGIN {{OFS=\"\\t\"}} $1 == \"{params.sg}\" {{ print $2,$3,$4}}' {input.tsv} > {output.tsv} 2> {log}; "
		" if [ $(wc -l {output.tsv} | awk '{{print $1}}') = \"0\" ]; then "
		"    echo 'No lines in regions file. Check that requested scaff group exists' >> {log}; exit 1; fi "
