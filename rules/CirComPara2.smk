rule generate_meta_vars_circompara2:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        meta_csv = "results/{sample}/CirComPara2/meta.csv",
        vars_py = "results/{sample}/CirComPara2/vars.py"
    params:
        sample_name = "{sample}",
        meta = "meta.csv",
        fa = "ref/reference/hg19.fa",
        gtf = "ref/annotation/hg19.gtf",
        segemehl = "ref/circompara2/segemehl/hg19.idx",
        bwa = "ref/circompara2/bwa/hg19",
        bowtie2 = "ref/circompara2/bowtie2/hg19",
        bowtie = "ref/circompara2/bowtie/hg19",
        star = "ref/circompara2/star/hg19",
        hisat2 = "ref/circompara2/hisat2/hg19",
        genepred = "ref/circompara2/hg19.genePred.wgn"
    script:
        os.path.join(workflow.basedir, "scripts/CirComPara2/generate_config.py")


rule detect_circompara2:
    input:
        meta_csv = "results/{sample}/CirComPara2/meta.csv",
        vars_py = "results/{sample}/CirComPara2/vars.py"
    output:
        #stats = "results/{sample}/CirComPara2/read_statistics/read_stats_collect/read_stats_collect.txt",
        time = "results/{sample}/CirComPara2/time.txt"
    container:
        "container/circompara2_new.sif"
    params:
        outdir = "results/{sample}/CirComPara2",
        time = "time.txt"
    resources:
        mem_mb = 64000,
        time = 720
    shell:
        """
        cd {params.outdir}
        /usr/bin/time -v -o {params.time} /circompara2/circompara2
        cd ../../..
        """