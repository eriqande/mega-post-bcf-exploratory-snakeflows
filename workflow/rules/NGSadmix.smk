


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
		" (CURDIR=$PWD && TMP=$(mktemp -d) && echo $TMP && cd $TMP && "
		" git clone {params.url} && "
        " cd NGSadmix  && "
        " git checkout {params.hash} && "
		" g++ NGSadmix.cpp -O3 -lpthread -lz -o $(basename {output}) && "
		" mv $(basename {output}) $CURDIR/{output} )  > {log} 2>&1 "





rule run_ngsadmix:
	input:
		beagle="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/beagle-gl.gz",
		bin="results/bin/NGSadmix"
	output:
		multiext("results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/ngsadmix/maf_{min_maf}/K_{K}_rep_{rep}/output", ".filter", ".fopt.gz", ".log", ".qopt")
	threads: 4
	resources:
		mem_mb=19200
	log:
		"results/logs/run_ngsadmix/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/maf_{min_maf}/K_{K}_rep_{rep}/log.txt"
	shell:
		" (OUTDIR={output[0]} && OUTDIR=${{OUTDIR/.filter/}} && "
		" {input.bin} -likes {input.beagle} -K {wildcards.K} -o $OUTDIR -P {threads} -minMaf {wildcards.min_maf}) > {log} 2>&1 "