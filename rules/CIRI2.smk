rule mapping_CIRI2:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        sam = "results/{sample}.{replicate}/CIRI2/{sample}.sam",
        log = "results/{sample}.{replicate}/CIRI2/{sample}_mapping.log",
        time = "results/{sample}.{replicate}/CIRI2/time1.txt"
    params:
        bwa = "ref/bwa/hg19",
        time = "results/{sample}.{replicate}/CIRI2/time1.txt"
    conda:
        os.path.join(workflow.basedir, "envs/CIRI2.yaml")
    threads:
        16
    resources:
        mem_mb = 48000,
        runtime = 720
    shell:
        "/usr/bin/time -v -o {params.time} bwa mem -T {threads} {params.bwa} {input.reads1} {input.reads2} 1> {output.sam} 2> {output.log}"

rule detect_CIRI2:
    input:
        sam = "results/{sample}.{replicate}/CIRI2/{sample}.sam"
    output:
        outfile = "results/{sample}.{replicate}/CIRI2/{sample}.out",
        log = "results/{sample}.{replicate}/CIRI2/{sample}.out.log",
        time = "results/{sample}.{replicate}/CIRI2/time2.txt"
    params:
        fa = "ref/reference/hg19.fa",
        gtf = "ref/annotation/hg19.gtf",
        time = "results/{sample}.{replicate}/CIRI2/time2.txt"
    conda:
        os.path.join(workflow.basedir, "envs/CIRI2.yaml")
    resources:
        mem_mb = 48000,
        runtime = 720
    shell:
        "/usr/bin/time -v -o {params.time} perl scripts/CIRI2/CIRI2.pl -I {input.sam} -O {output.outfile} -F {params.fa} -A {params.gtf}"

rule compile_time_CIRI2:
    input:
        time1 = "results/{sample}.{replicate}/CIRI2/time1.txt",
        time2 = "results/{sample}.{replicate}/CIRI2/time2.txt",
    output:
        time = "results/{sample}.{replicate}/CIRI2/time.txt"
    script:
        os.path.join(workflow.basedir, "scripts/time_parser/time_compiler.py")