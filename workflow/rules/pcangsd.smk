



# this is a rule that will build and install pcangsd into
# the active conda env.
rule install_pcangsd:
	params:
		hash=config["pcangsd_version"],
		url=config["pcangsd_url"]
	output:  
		flagfile=touch("results/flags/pcangsd_installed")
	conda:
		"../envs/pcangsd.yaml"
	log:
		"results/logs/install_pcangsd/log.txt"
	shell:
		"(TMP=$(mktemp -d) && cd $TMP && "
		" git clone {params.url} && "
        " cd pcangsd  && "
        " git checkout {params.hash} && "
        " python setup.py build_ext --inplace && "
        " pip3 install -e .  ) > {log} 2>&1  "




# then, whenever you need to use pcangsd, you call the same
# conda environment, and have the flagfile as an input
# depenency to ensure that pcangsd has already been successfully
# built into that conda env.
# the active conda env.  Yep! That works nicely.


# this one spits out the genotype posteriors and then beagle-izes them
rule pcangsd_with_gposts:
	input:  
		flagfile="results/flags/pcangsd_installed",
		beagle="results/beagle-gl/{bcf_id}/thin_{thin_int}_{thin_start}/beagle-gl.gz"
	output:
		cov="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/out.cov",
		gposts="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/out.gposts.tsv"
	conda:
		"../envs/pcangsd.yaml"
	threads: 20
	log:
		"results/logs/pcangsd_with_gposts/bcf_{bcf_id}/thin_{thin_int}_{thin_start}/log.txt"
	shell:
		"OUTPRE=$(dirname {output.gposts})/out && "
		"pcangsd -b {input.beagle} --minMaf 0.0 -t {threads} --post_save --out $OUTPRE > {log} 2>&1 "




