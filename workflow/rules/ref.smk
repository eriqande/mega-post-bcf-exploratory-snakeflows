



# this rule simply takes the fai from the fai_path
# in the config, and it puts a new one in the same
# location with the suffix -ANGSD.  This one has
# chromosome names in which the undescores have been
# replaced with dashes.  When this rule is run, it additionlly
# creates two files: prefix-angsd-names-second.tsh and prefix-angsd-names-first.tsv
# that would be suitable for changing the chromosome names in a bcf
# file using bcftools rename-chrs.
rule angsd_chromosome_names:
    input:
        config["fai_path"]
    output:
        fai="{p}-ANGSD".format(p=config["fai_path"]),
        angsd_second="{p}-angsd-names-second.tsv".format(p=config["fai_path"]),
        angsd_first="{p}-angsd-names-first.tsv".format(p=config["fai_path"]),
    log:
        "results/logs/angsd_chromosome_names/log.txt"
    shell:
        " (                                      "
        " awk '                                  "
        "   BEGIN {{OFS=\"\\t\"}}                "
        "   {{ gsub(/_/, \"-\", $1); print $0}}  "
        " ' {input} > {output.fai} &&            "
        "  paste {input} {output.fai} | cut -f 1,6 > {output.angsd_second} && "
        "  paste {output.fai} {input} | cut -f 1,6 > {output.angsd_first}     "
        " )  2> {log}                            "


