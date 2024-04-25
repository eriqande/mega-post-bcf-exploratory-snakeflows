



rule bcf_thin_scatter:
    input:
        bcf="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_0_0/main.bcf",
        pos="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_0_0/sections/{scaff_grp}.positions.tsv.gz"
    params:
        thin_int="{thin_int}",
        thin_start="{thin_start}",
        sg="{scaff_grp}"
    output:
        bcf="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.bcf",
        stats="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.bcf_stats.txt",
        pos="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.positions.tsv",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/bcf_thin_scatter/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.log"
    benchmark:
        "results/benchmarks/bcf_thin_scatter/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/sections/{scaff_grp}.bmk"
    shell:
        " (          "
        "    gunzip -c {input.pos} | awk -v thin_int={params.thin_int} -v thin_start={params.thin_start} -f workflow/scripts/thin_positions.awk > {output.pos} && "
        "    bcftools view -Ob -R {output.pos} {input.bcf}  > {output.bcf}  && "
        "    bcftools stats --af-tag MAF {output.bcf} > {output.stats} "
        " ) 2> {log} "




rule bcf_thin_gather:
    input:
        bcfs=expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_{{thin_int}}_{{thin_start}}/sections/{sg}.bcf", sg=unique_scaff_groups),
        poses=expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_{{thin_int}}_{{thin_start}}/sections/{sg}.positions.tsv", sg=unique_scaff_groups),
        statses=expand("results/bcf_{{bcf_id}}/filt_{{bcfilt}}/{{sampsub}}/thin_{{thin_int}}_{{thin_start}}/sections/{sg}.bcf_stats.txt", sg=unique_scaff_groups),
    output:
        bcf="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/main.bcf",
        csi="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/main.bcf.csi",
        pos="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/info/positions.tsv.gz",
        stats="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/info/bcf_stats.txt",
        samps="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/info/samples.txt",
    conda:
        "../envs/bcftools.yaml"
    log:
        "results/logs/bcf_thin_gather/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/main.log"
    benchmark:
        "results/benchmarks/bcf_thin_gather/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/main.bmk"
    shell:
        "( bcftools concat --naive {input.bcfs} > {output.bcf}     && "
        "  bcftools index {output.bcf}                             && "
        "  cat {input.poses} | gzip -c  > {output.pos}             && "
        "  plot-vcfstats -m {input.statses} > {output.stats}       && "
        "  bcftools query -l {output.bcf} > {output.samps}            "
        ") 2> {log} "

