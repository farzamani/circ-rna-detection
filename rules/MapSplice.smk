rule fastq_to_fasta_mapsplice:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        reads1 = temp("results/{sample}.{replicate}/MapSplice/{sample}_1.fa"),
        reads2 = temp("results/{sample}.{replicate}/MapSplice/{sample}_2.fa")
    resources:
        mem_mb = 4000,
        runtime = 720
    shell:
        """
        seqtk seq -a {input.reads1} > {output.reads1}
        seqtk seq -a {input.reads2} > {output.reads2}
        """

rule detect_mapsplice:
    input:
        reads1 = "results/{sample}.{replicate}/MapSplice/{sample}_1.fa",
        reads2 = "results/{sample}.{replicate}/MapSplice/{sample}_2.fa"
    output:
        outdir = directory("results/{sample}.{replicate}/MapSplice/mps3_results/"),
        time = "results/{sample}.{replicate}/MapSplice/time.txt"
    threads:
        16
    params:
        index = "ref/chrom_index",
        time = "results/{sample}.{replicate}/MapSplice/time.txt"
    resources:
        mem_mb = 72000,
        runtime = 1000
    shell:
        "/usr/bin/time -v -o {params.time} mps_regular_circRNA -G {params.index} -1 {input.reads1} -2 {input.reads2} -T {threads} -O {output.outdir}"