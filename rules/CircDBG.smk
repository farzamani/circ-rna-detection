rule config_script_circdbg:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        config = "results/{sample}.{replicate}/CircDBG/config.ini"
    params:
        gtf = "ref/annotation/hg19.gtf",
        fa = "ref/reference/hg19.fa"
    script:
        os.path.join(workflow.basedir, "scripts/CircDBG/circdbg.py")

rule detect_circdbg:
    input:
        config = "results/{sample}.{replicate}/CircDBG/config.ini"
    output:
        result = "results/{sample}.{replicate}/CircDBG/Detection_Result/Brief_sum.txt",
        time = "results/{sample}.{replicate}/CircDBG/time.txt"
    params:
        outdir = "results/{sample}.{replicate}/CircDBG",
        time = "time.txt"
    threads:
        16
    resources:
        mem_mb = 84000,
        runtime = 720
    shell:
        """
        cd {params.outdir}
        /usr/bin/time -v -o {params.time} CircRNADBG ./config.ini
        cd ../../..
        """