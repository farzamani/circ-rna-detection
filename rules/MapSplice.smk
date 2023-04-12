rule fastq_to_fasta_mapsplice:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        reads1 = temp("results/{sample}/MapSplice/{sample}_1.fa"),
        reads2 = temp("results/{sample}/MapSplice/{sample}_2.fa")
    resources:
        mem_mb = 4000,
        time = 720
    shell:
        """
        seqtk seq -a {input.reads1} > {output.reads1}
        seqtk seq -a {input.reads2} > {output.reads2}
        """

rule detect_mapsplice:
    input:
        reads1 = "results/{sample}/MapSplice/{sample}_1.fa",
        reads2 = "results/{sample}/MapSplice/{sample}_2.fa"
    output:
        outdir = directory("results/{sample}/MapSplice/mps3_results/"),
        time = "results/{sample}/MapSplice/time.txt"
    params:
        index = "ref/chrom_index",
        thread = 16,
        time = "results/{sample}/MapSplice/time.txt"
    resources:
        mem_mb = 12000,
        time = 720
    shell:
        "/usr/bin/time -v -o {params.time} mps_regular_circRNA -G {params.index} -1 {input.reads1} -2 {input.reads2} -T {params.thread} -O {output.outdir}"