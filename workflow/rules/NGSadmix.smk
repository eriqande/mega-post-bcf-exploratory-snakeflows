


# this is a rule that will build ngsadmix and put it somewhere simple for use
rule install_ngsadmix:
	params:
		hash=config["ngsAdmix"]["version"],
		url=config["ngsAdmix"]["url"]
	output:  
		"results/bin/NGSadmix"
	log:
		"results/logs/install_ngsadmix/log.txt"
	shell:
		"(TMP=$(mktemp -d) && echo $TMP && cd $TMP && "
		" git clone {params.url} && "
        " cd NGSadmix  && "
        " git checkout {params.hash} && "
		" g++ NGSadmix.cpp -O3 -lpthread -lz -o {output} "





rule run_ngsadmix:
	input:
		beagle="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/beagle-gl.gz",
		bin="results/bin/NGSadmix"
	output:
		"results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/ngsadmix/K_{K}_rep_{rep}/"
	threads: 4
	resources:
		mem_mb=19200
	shell:
		" ./NGSadmix -likes {input} -K {wildcards.K} -o {output} -P {threads} "