rule mapping_mate1_circtools:
    input:
        reads1 = rules.trim.output.trimmed1
    output:
        mate1 = directory("results/{sample}/DCC/mate1//"),
        junction1 = "results/{sample}/DCC/mate1/Chimeric.out.junction",
        time = "results/{sample}/DCC/time1.txt"
    params:
        star = "ref/dcc",
        gtf = "ref/annotation/hg19.gtf",
        outdir = directory("results/{sample}/DCC/mate1//"),
        time = "results/{sample}/DCC/time1.txt"
    conda:
        os.path.join(workflow.basedir, "envs/dcc.yaml")
    resources:
        mem_mb = 42000,
        time = 720
    shell:
        """
        /usr/bin/time -v -o {params.time} \
        STAR    --runThreadN 10\
                --genomeDir {params.star}\
                --genomeLoad NoSharedMemory\
                --readFilesIn {input.reads1}\
                --outFileNamePrefix {params.outdir}\
                --readFilesCommand zcat\
                --outReadsUnmapped Fastx\
                --outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
                --outSJfilterOverhangMin 15   15   15   15\
                --outFilterMultimapNmax 20\
                --outFilterScoreMin 1\
                --outFilterMatchNminOverLread 0.7\
                --outFilterMismatchNmax 999\
                --outFilterMismatchNoverLmax 0.05\
                --alignIntronMin 20\
                --alignIntronMax 1000000\
                --alignMatesGapMax 1000000\
                --alignSJoverhangMin 15\
                --alignSJDBoverhangMin 10\
                --alignSoftClipAtReferenceEnds No\
                --chimSegmentMin 15\
                --chimScoreMin 15\
                --chimScoreSeparation 10\
                --chimJunctionOverhangMin 15\
                --sjdbGTFfile {params.gtf}\
                --quantMode GeneCounts\
                --twopassMode Basic\
                --chimOutType Junctions SeparateSAMold
        """

rule mapping_mate2_circtools:
    input:
        reads2 = rules.trim.output.trimmed2
    output:
        mate2 = directory("results/{sample}/DCC/mate2//"),
        junction2 = "results/{sample}/DCC/mate2/Chimeric.out.junction",
        time = "results/{sample}/DCC/time2.txt"
    params:
        star = "ref/dcc",
        gtf = "ref/annotation/hg19.gtf",
        outdir = directory("results/{sample}/DCC/mate2//"),
        time = "results/{sample}/DCC/time2.txt"
    conda:
        os.path.join(workflow.basedir, "envs/dcc.yaml")
    resources:
        mem_mb = 42000,
        time = 720
    shell:
        """
        /usr/bin/time -v -o {params.time} \
        STAR    --runThreadN 10\
                --genomeDir {params.star}\
                --genomeLoad NoSharedMemory\
                --readFilesIn {input.reads2}\
                --outFileNamePrefix {params.outdir}\
                --readFilesCommand zcat\
                --outReadsUnmapped Fastx\
                --outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
                --outSJfilterOverhangMin 15   15   15   15\
                --outFilterMultimapNmax 20\
                --outFilterScoreMin 1\
                --outFilterMatchNminOverLread 0.7\
                --outFilterMismatchNmax 999\
                --outFilterMismatchNoverLmax 0.05\
                --alignIntronMin 20\
                --alignIntronMax 1000000\
                --alignMatesGapMax 1000000\
                --alignSJoverhangMin 15\
                --alignSJDBoverhangMin 10\
                --alignSoftClipAtReferenceEnds No\
                --chimSegmentMin 15\
                --chimScoreMin 15\
                --chimScoreSeparation 10\
                --chimJunctionOverhangMin 15\
                --sjdbGTFfile {params.gtf}\
                --quantMode GeneCounts\
                --twopassMode Basic\
                --chimOutType Junctions SeparateSAMold
        """

rule mapping_both_circtools:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        junction = "results/{sample}/DCC/samplesheet/Chimeric.out.junction",
        samplesheet = directory("results/{sample}/DCC/samplesheet//"),
        time = "results/{sample}/DCC/time3.txt"
    params:
        star = "ref/dcc",
        gtf = "ref/annotation/hg19.gtf",
        outdir = directory("results/{sample}/DCC/samplesheet//"),
        time = "results/{sample}/DCC/time3.txt"
    conda:
        os.path.join(workflow.basedir, "envs/dcc.yaml")
    resources:
        mem_mb = 42000,
        time = 720
    shell:
        """
        /usr/bin/time -v -o {params.time} \
        STAR    --runThreadN 10\
                --genomeDir {params.star}\
                --genomeLoad NoSharedMemory\
                --readFilesIn {input.reads1} {input.reads2}\
                --outFileNamePrefix {params.outdir}\
                --readFilesCommand zcat\
                --outReadsUnmapped Fastx\
                --outSAMattributes NH   HI   AS   nM   NM   MD   jM   jI   XS\
                --outSJfilterOverhangMin 15   15   15   15\
                --outFilterMultimapNmax 20\
                --outFilterScoreMin 1\
                --outFilterMatchNminOverLread 0.7\
                --outFilterMismatchNmax 999\
                --outFilterMismatchNoverLmax 0.05\
                --alignIntronMin 20\
                --alignIntronMax 1000000\
                --alignMatesGapMax 1000000\
                --alignSJoverhangMin 15\
                --alignSJDBoverhangMin 10\
                --alignSoftClipAtReferenceEnds No\
                --chimSegmentMin 15\
                --chimScoreMin 15\
                --chimScoreSeparation 10\
                --chimJunctionOverhangMin 15\
                --sjdbGTFfile {params.gtf}\
                --quantMode GeneCounts\
                --twopassMode Basic\
                --chimOutType Junctions SeparateSAMold
        """

rule detect_circtools:
    input:
        mate1 = "results/{sample}/DCC/mate1/Chimeric.out.junction",
        mate2 = "results/{sample}/DCC/mate2/Chimeric.out.junction",
        samplesheet = "results/{sample}/DCC/samplesheet/Chimeric.out.junction"
    output:
        circ = "results/{sample}/DCC/detect/CircCoordinates",
        rna_count = "results/{sample}/DCC/detect/CircRNACount",
        time = "results/{sample}/DCC/time4.txt"
    params:
        gtf = "ref/annotation/hg19.gtf",
        repeats = "ref/annotation/repeats.gtf",
        fa = "ref/reference/hg19.fa",
        outdir = directory("results/{sample}/DCC/detect//"),
        tmp = directory("results/{sample}/DCC/detect/temp//"),
        time = "results/{sample}/DCC/time4.txt"
    conda:
        os.path.join(workflow.basedir, "envs/dcc.yaml")
    resources:
        mem_mb = 42000,
        time = 720
    shell:
        """
        /usr/bin/time -v -o {params.time} \
        DCC {input.samplesheet} \
            -mt1 {input.mate1} \
            -mt2 {input.mate2} \
            -D \
            -T 20 \
            -R {params.repeats} \
            -an {params.gtf} \
            -Pi \
            -M \
            -fg \
            -A {params.fa} \
            -O {params.outdir} \
            -t {params.tmp}
        
        rm -rf {params.tmp}
        """

rule compile_time_circtools:
    input:
        time1 = "results/{sample}/DCC/time1.txt",
        time2 = "results/{sample}/DCC/time2.txt",
        time3 = "results/{sample}/DCC/time3.txt",
        time4 = "results/{sample}/DCC/time4.txt",
    output:
        time = "results/{sample}/DCC/time.txt"
    script:
        os.path.join(workflow.basedir, "scripts/time_parser/time_compiler.py")