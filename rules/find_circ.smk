rule make_bam_find_circ:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        bam = temp("results/{sample}.{replicate}/find_circ/{sample}.bam"),
        bowtie_log = "results/{sample}.{replicate}/find_circ/{sample}_bowtie.log",
        time = "results/{sample}.{replicate}/find_circ/time1.txt"
    params:
        bowtie2 = "ref/bowtie2/hg19",
        time = "results/{sample}.{replicate}/find_circ/time1.txt"
    conda:
        os.path.join(workflow.basedir, "envs/find_circ.yaml")
    threads:
        16
    resources:
        mem_mb = 12000,
        runtime = 720
    shell:
        "/usr/bin/time -v -o {params.time} bash -c 'bowtie2 -p {threads} --very-sensitive --score-min=C,-15,0 --mm -x {params.bowtie2} -q -1 {input.reads1} -2 {input.reads1} 2> {output.bowtie_log} | samtools view -hbuS - | samtools sort - > {output.bam}'"

rule unmapped_find_circ:
    input:
        bam = "results/{sample}.{replicate}/find_circ/{sample}.bam"
    output:
        unmapped = temp("results/{sample}.{replicate}/find_circ/{sample}_unmapped.bam"),
        time = "results/{sample}.{replicate}/find_circ/time2.txt"
    conda:
        os.path.join(workflow.basedir, "envs/find_circ.yaml")
    params:
        time = "results/{sample}.{replicate}/find_circ/time2.txt"
    resources:
        mem_mb = 12000,
        runtime = 720
    shell:
        "/usr/bin/time -v -o {params.time} bash -c 'samtools view -hf 4 {input.bam} | samtools view -Sb - > {output.unmapped}'"
        
rule unmapped2anchors_find_circ:
    input:
        unmapped = "results/{sample}.{replicate}/find_circ/{sample}_unmapped.bam" 
    output:
        fq = temp("results/{sample}.{replicate}/find_circ/{sample}.fastq.gz"),
        time = "results/{sample}.{replicate}/find_circ/time3.txt"
    conda:
        os.path.join(workflow.basedir, "envs/find_circ.yaml")
    params:
        time = "results/{sample}.{replicate}/find_circ/time3.txt"
    resources:
        mem_mb = 12000,
        runtime = 720
    shell:
        "/usr/bin/time -v -o {params.time} bash -c 'python2 scripts/find_circ/unmapped2anchors.py {input.unmapped} | gzip > {output.fq}'"

rule detect_find_circ:
    input:
        fq = "results/{sample}.{replicate}/find_circ/{sample}.fastq.gz",
    output:
        stats = "results/{sample}.{replicate}/find_circ/stats.log",
        fa = "results/{sample}.{replicate}/find_circ/spliced_reads.fa",
        bed = "results/{sample}.{replicate}/find_circ/spliced_sites.bed",
        time = "results/{sample}.{replicate}/find_circ/time4.txt"
    params:
        bowtie2 = "ref/bowtie2/hg19",
        ref = "ref/reference/hg19.fa",
        time = "results/{sample}.{replicate}/find_circ/time4.txt"
    conda:
        os.path.join(workflow.basedir, "envs/find_circ.yaml")
    threads:
        16
    resources:
        mem_mb = 48000,
        runtime = 720
    shell:
        "/usr/bin/time -v -o {params.time} bash -c 'bowtie2 -p {threads} --score-min=C,-15,0 --reorder --mm -q -U {input.fq} -x {params.bowtie2} | python2 scripts/find_circ/find_circ.py --genome={params.ref} --stats={output.stats} --reads={output.fa} > {output.bed}'"

# original shell: grep CIRCULAR {input.bed} | grep -v chrM | awk '$5>=2' | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | python2 scripts/find_circ/maxlength.py 100000 > {output.circ}
rule filter_find_circ:
    input:
        bed = "results/{sample}.{replicate}/find_circ/spliced_sites.bed"
    output:
        circ = "results/{sample}.{replicate}/find_circ/circ_candidates.bed",
        time = "results/{sample}.{replicate}/find_circ/time5.txt"
    conda:
        os.path.join(workflow.basedir, "envs/find_circ.yaml")
    params:
        time = "results/{sample}.{replicate}/find_circ/time5.txt"
    resources:
        mem_mb = 24000,
        runtime = 720
    shell:
        "/usr/bin/time -v -o {params.time} bash -c 'grep CIRCULAR {input.bed} | grep -v chrM | awk -f scripts/find_circ/filter.awk | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | python2 scripts/find_circ/maxlength.py 100000 > {output.circ}'"

rule compile_time_find_circ:
    input:
        time1 = "results/{sample}.{replicate}/find_circ/time1.txt",
        time2 = "results/{sample}.{replicate}/find_circ/time2.txt",
        time3 = "results/{sample}.{replicate}/find_circ/time3.txt",
        time4 = "results/{sample}.{replicate}/find_circ/time4.txt",
        time5 = "results/{sample}.{replicate}/find_circ/time5.txt"
    output:
        time = "results/{sample}.{replicate}/find_circ/time.txt"
    script:
        os.path.join(workflow.basedir, "scripts/time_parser/time_compiler.py")