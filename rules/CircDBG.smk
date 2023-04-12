rule config_script_circdbg:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        config = "results/{sample}/CircDBG/config.ini"
    params:
        gtf = "ref/annotation/hg19.gtf",
        fa = "ref/reference/hg19.fa"
    script:
        os.path.join(workflow.basedir, "scripts/CircDBG/circdbg.py")

rule detect_circdbg:
    input:
        config = "results/{sample}/CircDBG/config.ini"
    output:
        result = "results/{sample}/CircDBG/Detection_Result/Brief_sum.txt",
        time = "results/{sample}/CircDBG/time.txt"
    params:
        outdir = "results/{sample}/CircDBG",
        time = "time.txt"
    resources:
        mem_mb = 12000,
        time = 720
    shell:
        """
        cd {params.outdir}
        /usr/bin/time -v -o {params.time} CircRNADBG ./config.ini
        cd ../../..
        """