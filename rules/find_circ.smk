def find_circ_v3_input(wildcards):
    output_dict = {}

    if config['filter']['UMI']:
        output_dict['pair1'] = rules.umi_collapse.output.pair1
        output_dict['pair2'] = rules.umi_collapse.output.pair2
    else:
        output_dict['pair1'] = rules.trim_adapter.params.pair1
        output_dict['pair2'] = rules.trim_adapter.params.pair2

    output_dict['BowtieIndexfile'] = config["find_circ_v3"]["BowtieIndexfile"]

    return output_dict

rule find_circ_v3:
    input:
        unpack(find_circ_v3_input)
    output:
        bamfile=temp("results/find_circ_output/{sample}.bam"),
        unmapped_bam=temp("results/find_circ_output/{sample}_unmapped.bam"),
        anchors=temp("results/find_circ_output/{sample}.anchors.qfa.gz"),
        sites_bed="results/find_circ_output/{sample}/{sample}.sample.sites.bed",
        sites_reads="results/find_circ_output/{sample}/{sample}.sites.reads",
        circ_candidates="results/find_circ_output/{sample}/{sample}.circ_candidates.bed",
        circ_candidates_35x2="results/find_circ_output/{sample}/{sample}.circ_candidates_map35x2.bed"
    params:
        BowtieIndex=config["find_circ"]["BowtieIndex"],
        genomedir=config["find_circ"]["genomedir"]
    log:
        bt2_firstpass="logs/{sample}.bt2_firstpass.log",
        sites="logs/find_circ_output/{sample}/{sample}.sites.log"
    conda:
        os.path.join(workflow.basedir, "envs/find_circ.yml")
    resources:
        mem_mb=12000,
        time=720
    threads:
        14
    shell:
        """
        echo "Running bowtie2 to find unmapped reads"
        ## Bowtie2 - find unmapped reads
        bowtie2 -p {threads} --very-sensitive --mm -M20 --score-min=C,-15,0 \
          -x {params.BowtieIndex} -q -1 {input.pair1} -2 {input.pair2} \
          2> {log.bt2_firstpass} | samtools view -hbuS - | samtools sort - > {output.bamfile}

        ## Extract unmapped reads with samtools
        echo "Extract unmapped reads"
        samtools view -hf 4 {output.bamfile} | samtools view -Sb - > {output.unmapped_bam}        

        echo "Producing anchors"
        python scripts/find_circ/unmapped2anchors.py {output.unmapped_bam} | gzip > {output.anchors}

        echo "Final bowtie"
        ### Bowtie2
        bowtie2 -p {threads} --reorder --mm -M20 --score-min=C,-15,0 -q -x {params.BowtieIndex} -U {output.anchors} |
          python scripts/find_circ/find_circ_v2.py -G {params.genomedir} -p {wildcards.sample} -s {log.sites} \
          > {output.sites_bed} 2> {output.sites_reads}

        echo "Applying thresholds"
        grep circ_ {output.sites_bed} | grep -v chrM |
            python scripts/find_circ/sum.py -2,3 |
            python scripts/find_circ/scorethresh.py -20 1 |
            python scripts/find_circ/scorethresh.py -19 2 |
            python scripts/find_circ/scorethresh.py -18 2 |
            python scripts/find_circ/scorethresh.py 7 2 |
            python scripts/find_circ/scorethresh.py 9,10 35 |
            python scripts/find_circ/scorethresh.py -21 10000 > {output.circ_candidates}
        
        echo "Applything thresholds with 35x2"
        grep circ_ {output.sites_bed} | grep -v chrM |
            python scripts/find_circ/sum.py -2,3 |
            python scripts/find_circ/scorethresh.py -20 1 |
            python scripts/find_circ/scorethresh.py -19 2 |
            python scripts/find_circ/scorethresh.py -18 2 |
            python scripts/find_circ/scorethresh.py 7 2 |
            python scripts/find_circ/scorethresh.py 9 35 |
            python scripts/find_circ/scorethresh.py 10 35 |
            python scripts/find_circ/scorethresh.py -21 10000 > {output.circ_candidates_35x2}

        """