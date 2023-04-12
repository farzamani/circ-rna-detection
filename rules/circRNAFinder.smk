rule convert_reads_circRNAFinder:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        reads = "results/{sample}/circRNAFinder/reads_{sample}.fa"
    conda:
        os.path.join(workflow.basedir, "envs/circRNAFinder.yaml")
    resources:
        mem_mb = 12000,
        time = 720
    shell:
        "scripts/circRNAFinder/fastx_collapser.py -f {input.reads1},{input.reads2} {output.reads}"

rule config_circRNAFinder:
    input:
        reads = "results/{sample}/circRNAFinder/reads_{sample}.fa"
    output:
        config_find = "results/{sample}/circRNAFinder/circRNAFind.cfg",
        config_anno = "results/{sample}/circRNAFinder/circRNAAnno.cfg"
    params:
        sample_name = "{sample}",
        gtf_build = "ref/annotation/hg19.gtf.build",
        fa = "ref/reference/hg19.fa",
        bowtie = "ref/bowtie/hg19",
        transcriptome = "ref/transcriptome/hg19_trans"
    script:
        os.path.join(workflow.basedir, "scripts/circRNAFinder/generate_config.py")

rule detect_circRNAFinder:
    input:
        config_find = "results/{sample}/circRNAFinder/circRNAFind.cfg"
    output:
        outdir = "results/{sample}/circRNAFinder/{sample}/output/{sample}.circ.txt",
        temfile = temp(directory("results/{sample}/circRNAFinder/{sample}/temp")),
        time = "results/{sample}/circRNAFinder/time.txt"
    conda:
        os.path.join(workflow.basedir, "envs/circRNAFinder.yaml")
    params:
        outdir = "results/{sample}/circRNAFinder",
        time = "results/{sample}/circRNAFinder/time.txt",
        sample_name = "{sample}"
    resources:
        mem_mb = 12000,
        time = 720
    shell:
        """
        export PATH=$PATH:scripts/circRNAFinder
        /usr/bin/time -v -o {params.time} circRNAFind.py {input.config_find} -s 0
        rsync -a {params.sample_name} {params.outdir}
        rm -rf {params.sample_name}
        """

        