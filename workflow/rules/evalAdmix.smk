rule install_evaladmix:
  params:
    hash=config["evalAdmix"]["version"], 
    url=config["evalAdmix"]["url"]
  output:
    "results/bin/evalAdmix"
  log:
    "results/logs/install_evaladmix/log.txt"
  shell:
    " (CURDIR=$PWD && TMP=$(mktemp -d) && echo $TMP && cd $TMP && "
    " git clone {params.url} && "
    " cd evalAdmix  && "
    " git checkout {params.hash} && "
    " g++ -c filereader_and_conversions.cpp -O3 && "
    " g++ -c extractors.cpp -O3 && "
    " gcc -c -fPIC alloc.cpp -O3 && "
    " gcc -c asort.cpp -O3 && "
    " gcc -c evalAdmix.cpp -O3 && "
    " gcc -c ngsevalAdmix.cpp -O3 && "
    " g++  -o evalAdmix Cinterface.cpp filereader_and_conversions.o extractors.o asort.o alloc.o evalAdmix.o ngsevalAdmix.o -O3 -lz -lpthread -o $(basename {output}) &&"
    " mv $(basename {output}) $CURDIR/{output} )  > {log} 2>&1 "
    

rule run_evaladmix:
  input:
    beagle="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/beagle-gl/beagle-gl.gz",
    fopt="results/bcf_{bcf_id}/filt_{bcfilt}/{subsamp}/thin_{thin_int}_{thin_start}/ngsadmix/maf_{min_maf}/K_{K}_rep_{rep}/output.fopt.gz", 
    qopt="results/bcf_{bcf_id}/filt_{bcfilt}/{subsamp}/thin_{thin_int}_{thin_start}/ngsadmix/maf_{min_maf}/K_{K}_rep_{rep}/output.qopt",
    bin="results/bin/evalAdmix"
  output:
    "results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/evaladmix/maf_{min_maf}/K_{K}_rep_{rep}/output.corres.txt"
  threads: 10
  resources:
    mem_mb=19200
  log:
    "results/logs/run_evaladmix/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/maf_{min_maf}/K_{K}_rep_{rep}/log.txt"
  shell:
    " {input.bin} -beagle {input.beagle} \ "
    " -fname {input.fopt} \ "
    " -qname {input.qopt} \ "
    " -P {threads} \ "
    " -o {output} > {log} 2>&1 "
    
    
# simple rule to add the sample names to the evalAdmix corres file
rule postprocess_evaladmix:
    input:
        corres="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/evaladmix/maf_{min_maf}/K_{K}_rep_{rep}/output.corres.txt",
        samp_list="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/info/samples.txt"
    output:
        corres="results/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/evaladmix/maf_{min_maf}/K_{K}_rep_{rep}/output.corres.txt_with_sample_names"
    params:
        K="{K}"
    threads: 1
    log:
        "results/logs/postprocess_evaladmix/bcf_{bcf_id}/filt_{bcfilt}/{sampsub}/thin_{thin_int}_{thin_start}/maf_{min_maf}/K_{K}_rep_{rep}/log.txt"
    shell:
        "paste {input.samp_list} {input.corres} | "
        " awk -v K={params.K} ' "
        "  BEGIN {{OFS=\"\\t\"; printf(\"sample\"); for(i=1;i<=K;i++) printf(\"\\tQ%d\", i); printf(\"\\n\")}} "
        " {{print}} ' > {output.corres} 2> {log} "
    