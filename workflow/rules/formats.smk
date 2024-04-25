# Rules related to switching between formats and things




# makes beagle GL file from the PL field in BCF file.  Note that we capitate them
# and keep them around in case the individual sections are useful down the road
rule bcf2beagle_gl_scatter:
    input:
        bcf="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/main.bcf",
        csi="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/main.bcf.csi",
        regions="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/scaff_members/{scaff_grp}.scaff_members.tsv",
        sfile="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/info/samples.txt"
    output:
        body=temp("results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/sections/{scaff_grp}.body.gz"),
        top_row=temp("results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/sections/{scaff_grp}.toprow.gz"),
        beag="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/sections/{scaff_grp}.beagle-gl.gz"
    log:
        "results/logs/bcf2beagle_gl_scatter/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/{scaff_grp}.log"
    benchmark:
        "results/benchmarks/bcf2beagle_gl_scatter/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/{scaff_grp}.bmk"
    conda:
        "../envs/bcftools.yaml"
    shell:
        " (           " 
        "    awk -f workflow/scripts/beagle3header.awk {input.sfile} | gzip -c > {output.top_row}  && "
        "    bcftools view -Ou -R {input.regions} {input.bcf} |  "
        "    bcftools query -f '%CHROM:%POS\\t%REF\\t%ALT[\\t%PL]\\n' | "
        "    awk -f workflow/scripts/pl2gl.awk | gzip -c  >  {output.body} && "
        "    cat {output.top_row} {output.body} > {output.beag}  " 
        " ) 2> {log}  "



rule bcf2beagle_gl_gather:
    input: 
        header=expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_{{thin_int}}_{{thin_start}}/beagle-gl/sections/{sg}.toprow.gz", sg=first_scaff_group_id),
        scaff_gzs=expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_{{thin_int}}_{{thin_start}}/beagle-gl/sections/{sg}.body.gz", sg=unique_scaff_groups)
    output:
        "results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/beagle-gl.gz"
    log:
        "results/logs/bcf2beagle_gl_gather/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/main.log"
    benchmark:
        "results/logs/bcf2beagle_gl_gather/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/main.bmk"
    shell:
        "cat {input.header} {input.scaff_gzs} > {output} 2> {log} "




