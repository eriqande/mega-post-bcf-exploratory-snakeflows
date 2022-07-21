



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
	params: 
		minMaf = "{min_maf}"
	output:
		args="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/out.args",
		cov="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/out.cov",
		gposts=temp("results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/out.gpost.tsv"),
		mafs="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/out.maf.npy",
		sites="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/out.sites"
		
	conda:
		"../envs/pcangsd.yaml"
	threads: 20
	log:
		pcangsd="results/logs/pcangsd_with_gposts/bcf_{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/pcangsd_part.txt",
		beagle="results/logs/pcangsd_with_gposts/bcf_{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/beagle_paste_part.txt",
	shell:
		" (OUTPRE=$(dirname {output.gposts})/out && "
		" pcangsd -b {input.beagle} --minMaf {params.minMaf} -t {threads} --post_save --maf_save --sites_save --out $OUTPRE > {log.pcangsd} 2>&1) "
		


rule pcangsd_beagle_post_bung:
	input:
		beagle="results/beagle-gl/{bcf_id}/thin_{thin_int}_{thin_start}/beagle-gl.gz",
		gposts="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/out.gpost.tsv",
		sites="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/out.sites"
	output:
		beagle_posts="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/beagle-post.gz",
		beagle_header="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/beagle_header"
	log:
		"results/logs/pcangsd_with_gposts/bcf_{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/pcangsd_beagle_post_bung/log.txt"
	shell:
		" (gunzip -c {input.beagle} | head -n 1 > {output.beagle_header} 2> {log})  && "
		" (gunzip -c {input.beagle}  | awk 'BEGIN {{OFS=\"\\t\"}} NR>1 {{print $1, $2, $3}}'  | "
		" paste {input.sites} - | awk 'BEGIN {{OFS=\"\\t\"}} $1==1 {{print $2, $3, $4}}' | "
		" paste - {input.gposts} | cat {output.beagle_header} - | gzip - >  {output.beagle_posts} 2>> {log}) "

