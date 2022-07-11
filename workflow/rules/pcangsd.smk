



# this is a rule that will build and install pcangsd into
# the active conda env.
rule install_pcangsd:
	params:
		hash=config["pcangsd_version"]
	output:  
		flagfile=touch("results/flags/pcangsd_installed")
	conda:
		"../envs/pcangsd.yaml"
	log:
		"results/logs/install_pcangsd/log.txt"
	shell:
		"(echo pcangsd version, git hash: {params.hash} > {log} && "
		" TMP=$(mktemp -d) && cd $TMP && "
		" git clone https://github.com/Rosemeis/pcangsd.git && "
        " cd pcangsd  && "
        " git checkout {params.hash} && "
        " python setup.py build_ext --inplace && "
        " pip3 install -e .  ) >> {log} 2>>&1  "

