



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
		"(TMP=$(mktemp -d) && cd $TMP && "
		" git clone https://github.com/Rosemeis/pcangsd.git && "
        " cd pcangsd  && "
        " git checkout {params.hash} && "
        " python setup.py build_ext --inplace && "
        " pip3 install -e .  ) > {log} 2>&1  "




# then, whenever you need to use pcangsd, you call the same
# conda environment, and have the flagfile as an input
# depenency to ensure that pcangsd has already been successfully
# built into that conda env.
# the active conda env.
rule print_pcangsd_message:
	input:  
		flagfile=touch("results/flags/pcangsd_installed")
	output:
		"Just-A-Test-File.txt"
	conda:
		"../envs/pcangsd.yaml"
	log:
		"results/logs/print_pcangsd_message/log.txt"
	shell:
		"pcangsd > {output} 2>&1 "