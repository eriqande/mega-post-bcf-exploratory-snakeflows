





# this is for simple pcangsd with no genotype posteriors.  It also computes the selection statistics.
rule pcangsd_no_gposts:
    input:  
        #flagfile="results/flags/pcangsd_installed",
        beagle="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/beagle-gl.gz"
    params: 
        minMaf = "{min_maf}"
    output:
        cov="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd_plain/maf_{min_maf}/{param_set}/out.cov",
        freqs="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd_plain/maf_{min_maf}/{param_set}/out.freqs",
        sites="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd_plain/maf_{min_maf}/{param_set}/out.sites"
        
    conda:
        "../envs/pcangsd.yaml"
    threads: 20
    resources:
        mem_mb=190000
    log:
        pcangsd="results/logs/pcangsd_no_gposts/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/maf_{min_maf}/{param_set}/pcangsd_part.txt",
        beagle="results/logs/pcangsd_no_gposts/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/maf_{min_maf}/{param_set}/beagle_paste_part.txt",
    shell:
        " (OUTPRE=$(dirname {output.cov})/out && "
        " pcangsd -b {input.beagle} --maf {params.minMaf} -t {threads} --maf-save --sites-save --selection --out $OUTPRE > {log.pcangsd} 2>&1) "



# this one spits out the genotype posteriors and then beagle-izes them
rule pcangsd_with_gposts:
    input:  
        #flagfile="results/flags/pcangsd_installed",
        beagle="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/beagle-gl.gz"
    params: 
        minMaf = "{min_maf}"
    output:
        cov="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/out.cov",
        gposts="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/out.post",
        freqs="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/out.freqs",
        sites="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/out.sites"
        
    conda:
        "../envs/pcangsd.yaml"
    threads: 20
    resources:
        mem_mb=190000
    log:
        pcangsd="results/logs/pcangsd_with_gposts/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/maf_{min_maf}/pcangsd_part.txt",
        beagle="results/logs/pcangsd_with_gposts/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/maf_{min_maf}/beagle_paste_part.txt",
    shell:
        " (OUTPRE=$(dirname {output.gposts})/out && "
        " pcangsd -b {input.beagle} --maf {params.minMaf} -t {threads} --post --maf-save --sites-save --out $OUTPRE > {log.pcangsd} 2>&1) "
        


rule pcangsd_beagle_post_bung:
    input:
        beagle="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/beagle-gl.gz",
        gposts="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/out.post",
        sites="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/out.sites",
    output:
        beagle_posts="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/beagle-post.gz",
        beagle_header="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/beagle_header"
    log:
        "results/logs/pcangsd_beagle_post_bung/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/maf_{min_maf}/log.txt"
    shell:
        " set +o pipefail; (gunzip -c {input.beagle} | awk 'NR==1 {{print; exit 0}}' > {output.beagle_header} 2> {log})  && "
        " (gunzip -c {input.beagle}  | awk 'BEGIN {{OFS=\"\\t\"}} NR>1 {{print $1, $2, $3}}'  | "
        " paste {input.sites} - | awk 'BEGIN {{OFS=\"\\t\"}} $1==1 {{print $2, $3, $4}}' | "
        " paste - {input.gposts} | cat {output.beagle_header} - | gzip - >  {output.beagle_posts} 2>> {log}) "




# this version of it takes the final beagle posterior output, and, instead of
# gzipping it, we run it through an awk script that breaks it out into
# a bunch of separate NON-GZIPPED beagle files, one for each scaffold group.
# The idea behind this is that we can then scatter angsd doAsso over scaffold groups.
rule pcangsd_beagle_post_slice:
    input:
        beagle="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/beagle-gl.gz",
        gposts="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/out.post",
        sites="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/out.sites",
        scaff_grp_path=get_scaff_group_path
    output:
        beagle_sections=expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_{{thin_int}}_{{thin_start}}/pcangsd/maf_{{min_maf}}/sections/{asg}-beagle-post.gz", asg = unique_scaff_groups),
        beagle_header="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/beagle_header"
    log:
        "results/logs/pcangsd_beagle_post_slice/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/maf_{min_maf}/log.txt"
    threads: 4
    shell:
        " set +o pipefail; (gunzip -c {input.beagle} | awk 'NR==1 {{print; exit 0}}' > {output.beagle_header} 2> {log})  && "
        " (gunzip -c {input.beagle}  | awk 'BEGIN {{OFS=\"\\t\"}} NR>1 {{print $1, $2, $3}}'  | "
        " paste {input.sites} - | awk 'BEGIN {{OFS=\"\\t\"}} $1==1 {{print $2, $3, $4}}' | "
        " paste - {input.gposts} | cat {input.scaff_grp_path} {output.beagle_header} - | "
        " awk -v path=\"$(dirname {output.beagle_header})/sections\" -v ext=post -f workflow/scripts/beagle-slicer.awk  >>  {log} 2>&1 ) "





# this is just a quick rule to extract sites from the beagle-post files
# This was for one-off use at one point...
rule extract_sites_from_beagle_posts:
    input:
        beag="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/sections/scaff_grp_{sgn}-beagle-post.gz",
        extr="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/extracto/scaff_grp_{sgn}.txt",
    output:
        extd="results/pcangsd/{bcf_id}/thin_{thin_int}_{thin_start}/maf_{min_maf}/extracted/scaff_grp_{sgn}-beag-posts.tsv",
    shell:
        " zcat {input.beag} | cat {input.extr} - | "
        " awk ' "
        "    BEGIN {{OFS=\"\t\"}} "
        "    NF==1 {{sg[$1]++; n++; next}} "
        "    /marker/ {{print; next}}  "
        "    $1 in sg {{print; m++; if(m==n) exit 0}} ' > {output.extd} "




# this creates a dot-samples file that has the PCs in it
rule attach_PCs_to_dotsamples:
    input:
        dots=lambda wc: config["bcf"][wc.bcf_id]["sample_subsets"][wc.sampsub]["dotsample"],
        cov="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/out.cov"
    params:
        num_pcs="{npc}"
    log:
        "results/logs/attach_PCs_to_dotsamples/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/npc_{npc}.log"
    envmodules:
        config["Rmodule"]
    output:
        dots="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/pcangsd/maf_{min_maf}/dotsample_PCs_{npc}.tsv"
    script:
        "../scripts/attach-pcs-to-dotsamples.R"


















