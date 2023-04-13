rule mapping_circsplice:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        sam = "results/{sample}/CircSplice/star/Chimeric.out.sam",
        time = "results/{sample}/CircSplice/time1.txt"
    params:
        star = "ref/dcc",
        gtf = "ref/annotation/hg19.gtf",
        outdir = directory("results/{sample}/CircSplice/star//"),
        time = "results/{sample}/CircSplice/time1.txt"
    conda:
        os.path.join(workflow.basedir, "envs/CircSplice.yaml")
    resources:
        mem_mb = 42000,
        time = 720
    shell:
        """
        /usr/bin/time -v -o {params.time}\
        STAR    --genomeDir {params.star}\
                --readFilesIn {input.reads1} {input.reads2}\
                --sjdbGTFfile {params.gtf}\
                --readFilesCommand zcat\
                --runThreadN 10\
                --chimSegmentMin 20\
                --chimScoreMin 1\
                --alignIntronMax 100000\
                --outFilterMismatchNmax 4\
                --alignTranscriptsPerReadNmax 100000\
                --outFilterMultimapNmax 2\
                --outFileNamePrefix {params.outdir}\
                --chimOutType Junctions SeparateSAMold
        """

rule run_circsplice:
    input:
        sam = "results/{sample}/CircSplice/star/Chimeric.out.sam"
    output:
        result_as = "results/{sample}/CircSplice/star/Chimeric.out.sam.result.as",
        result_circ = "results/{sample}/CircSplice/star/Chimeric.out.sam.result.circ",
        time = "results/{sample}/CircSplice/time2.txt"
    params:
        gtf = "ref/annotation/hg19.gtf",
        refflat = "ref/bed-refFlat.txt",
        fa = "ref/reference/hg19.fa",
        time = "results/{sample}/CircSplice/time2.txt"
    conda:
        os.path.join(workflow.basedir, "envs/CircSplice.yaml")
    resources:
        mem_mb = 20000,
        time = 720
    shell:
        "/usr/bin/time -v -o {params.time} perl scripts/CircSplice/CircSplice.pl {input.sam} {params.fa} {params.refflat}"

rule compile_time_circsplice:
    input:
        time1 = "results/{sample}/CircSplice/time1.txt",
        time2 = "results/{sample}/CircSplice/time2.txt"
    output:
        time = "results/{sample}/CircSplice/time.txt"
    script:
        os.path.join(workflow.basedir, "scripts/time_parser/time_compiler.py")


