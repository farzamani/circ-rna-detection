rule config_circmarker:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        config = "results/{sample}.{replicate}/CircMarker/config.ini"
    params:
        gtf = "ref/annotation/hg19.gtf",
        fa = "ref/reference/hg19.fa",
        bwa = "ref/bwa/hg19",
        outdir = "results/{sample}.{replicate}/CircMarker"
    script:
        os.path.join(workflow.basedir, "scripts/CircMarker/circmarker.py")

rule detect_circmarker:
    input:
        config = "results/{sample}.{replicate}/CircMarker/config.ini"
    output:
        result = "results/{sample}.{replicate}/CircMarker/Detection_Result/Brief_sum.txt",
        time = "results/{sample}.{replicate}/CircMarker/time.txt"
    params:
        outdir = "results/{sample}.{replicate}/CircMarker",
        time = "time.txt"
    threads:
        16
    resources:
        mem_mb = 84000,
        runtime = 720
    shell:
        """
        cd {params.outdir}
        /usr/bin/time -v -o {params.time} CircMarker ./config.ini
        cd ../../..
        """