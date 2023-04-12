# First step
rule alignment_circexplorer2:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        outdir = temp(directory("results/{sample}/CIRCexplorer2/tophat_fusion")),
        bam = "results/{sample}/CIRCexplorer2/tophat_fusion/accepted_hits.bam",
        log = "results/{sample}/CIRCexplorer2/CIRCexplorer2_alignment.log",
        time = "results/{sample}/CIRCexplorer2/time1.txt"
    threads:
        15
    params:
        bowtie2 = "ref/bowtie2/hg19",
        time = "results/{sample}/CIRCexplorer2/time1.txt"
    conda:
        os.path.join(workflow.basedir, "envs/circexplorer2.yaml") 
    resources:
        mem_mb = 12000,
        time = 720
    shell:
        "/usr/bin/time -v -o {params.time} tophat2 -o {output.outdir} -p {threads} --fusion-search --keep-fasta-order --no-coverage-search {params.bowtie2} {input.reads1} {input.reads1} > {output.log}"

# Second first-step
rule parse_circexplorer2:
    input:
        bam = "results/{sample}/CIRCexplorer2/tophat_fusion/accepted_hits.bam"
    output:
        bed = "results/{sample}/CIRCexplorer2/back_spliced_junction.bed",
        log = "results/{sample}/CIRCexplorer2/CIRCexplorer2_parse.log",
        time = "results/{sample}/CIRCexplorer2/time2.txt"
    params:
        aligner = "TopHat-Fusion",
        time = "results/{sample}/CIRCexplorer2/time2.txt"
    conda:
        os.path.join(workflow.basedir, "envs/circexplorer2.yaml")
    resources:
        mem_mb = 12000,
        time = 720
    shell:
        "/usr/bin/time -v -o {params.time} CIRCexplorer2 parse --pe -t {params.aligner} {input.bam} -b {output.bed} > {output.log}"

# Second step
rule annotate_circexplorer2:
    input:
        r = "ref/annotation/hg19.txt",
        g = "ref/reference/hg19.fa",
        bed = "results/{sample}/CIRCexplorer2/back_spliced_junction.bed"
    output:
        circ_known = "results/{sample}/CIRCexplorer2/circularRNA_known.txt",
        log = "results/{sample}/CIRCexplorer2/CIRCexplorer2_annotate.log",
        time = "results/{sample}/CIRCexplorer2/time3.txt"
    conda:
        os.path.join(workflow.basedir, "envs/circexplorer2.yaml")
    params:
        time = "results/{sample}/CIRCexplorer2/time3.txt"
    resources:
        mem_mb = 12000,
        time = 720
    shell:
        "/usr/bin/time -v -o {params.time} CIRCexplorer2 annotate -r {input.r} -g {input.g} -b {input.bed} -o {output.circ_known} > {output.log}"

# Third step
rule assembly_circexplorer2:
    input:
        r = "ref/annotation/hg19.txt",
        tophat_dir = "results/{sample}/CIRCexplorer2/tophat_fusion"
    output:
        assemble = temp(directory("results/{sample}/CIRCexplorer2/assemble")),
        log = "results/{sample}/CIRCexplorer2/CIRCexplorer2_assemble.log",
        time = "results/{sample}/CIRCexplorer2/time4.txt"
    conda:
        os.path.join(workflow.basedir, "envs/circexplorer2.yaml")
    params:
        time = "results/{sample}/CIRCexplorer2/time4.txt"
    resources:
        mem_mb = 12000,
        time = 720
    shell:
        "/usr/bin/time -v -o {params.time} CIRCexplorer2 assemble -r {input.r} -m {input.tophat_dir} -o {output.assemble} > {output.log}"

# Fourth step
rule denovo_circexplorer2:
    input:
        r = "ref/annotation/hg19.txt",
        g = "ref/reference/hg19.fa",
        tophat_dir = "results/{sample}/CIRCexplorer2/tophat_fusion",
        bed = "results/{sample}/CIRCexplorer2/back_spliced_junction.bed",
        assemble = "results/{sample}/CIRCexplorer2/assemble"
    output:
        log = "results/{sample}/CIRCexplorer2/CIRCexplorer2_denovo.log",
        abs = directory("results/{sample}/CIRCexplorer2/abs"),
        denovo = directory("results/{sample}/CIRCexplorer2/denovo"),
        time = "results/{sample}/CIRCexplorer2/time5.txt"
    conda:
        os.path.join(workflow.basedir, "envs/circexplorer2.yaml")
    params:
        time = "results/{sample}/CIRCexplorer2/time5.txt"
    resources:
        mem_mb = 12000,
        time = 720
    shell:
        "/usr/bin/time -v -o {params.time} CIRCexplorer2 denovo --abs {output.abs} -m {input.tophat_dir} -r {input.r} -g {input.g} -b {input.bed} -d {input.assemble} -o {output.denovo} > {output.log}"

# Compile time
rule compile_time_circexplorer2:
    input:
        time1 = "results/{sample}/CIRCexplorer2/time1.txt",
        time2 = "results/{sample}/CIRCexplorer2/time2.txt",
        time4 = "results/{sample}/CIRCexplorer2/time4.txt",
        time5 = "results/{sample}/CIRCexplorer2/time5.txt",
    output:
        time = "results/{sample}/CIRCexplorer2/time.txt"
    script:
        os.path.join(workflow.basedir, "scripts/time_parser/time_compiler.py")