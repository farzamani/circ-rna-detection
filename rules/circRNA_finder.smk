rule star_run_circRNA_finder:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2
    output:
        log = "results/{sample}.{replicate}/circRNA_finder/star_result/star_Log.out",
        junction = "results/{sample}.{replicate}/circRNA_finder/star_result/star_Chimeric.out.junction",
        outdir = directory("results/{sample}.{replicate}/circRNA_finder/star_result"),
        time = "results/{sample}.{replicate}/circRNA_finder/time1.txt"
    params:
        maxMismatch = 0.02,
        star = "ref/star",
        prefix = "results/{sample}.{replicate}/circRNA_finder/star_result/star_",
        time = "results/{sample}.{replicate}/circRNA_finder/time1.txt"
    conda:
        os.path.join(workflow.basedir, "envs/circRNA_finder.yaml")
    resources:
        mem_mb = 60000,
        runtime = 720
    shell:
        "/usr/bin/time -v -o {params.time} perl scripts/circRNA_finder/runStar.pl --inFile1 {input.reads1} --inFile2 {input.reads2} --genomeDir {params.star} --maxMismatch {params.maxMismatch} --outPrefix {params.prefix}"

rule detect_circRNA_finder:
    input:
        star = "results/{sample}.{replicate}/circRNA_finder/star_result//"
    output:
        junction = "results/{sample}.{replicate}/circRNA_finder/star_resultstar_filteredJunctions.bed",
        time = "results/{sample}.{replicate}/circRNA_finder/time2.txt"
    params:
        prefix = "results/{sample}.{replicate}/circRNA_finder/star_result",
        time = "results/{sample}.{replicate}/circRNA_finder/time2.txt"
    conda:
        os.path.join(workflow.basedir, "envs/circRNA_finder.yaml")
    resources:
        mem_mb = 60000,
        runtime = 720
    shell:
        "/usr/bin/time -v -o {params.time} perl scripts/circRNA_finder/postProcessStarAlignment.pl --starDir {input.star}/ --outDir {params.prefix}"
        
rule compile_time_circRNA_finder:
    input:
        time1 = "results/{sample}.{replicate}/circRNA_finder/time1.txt",
        time2 = "results/{sample}.{replicate}/circRNA_finder/time2.txt"
    output:
        time = "results/{sample}.{replicate}/circRNA_finder/time.txt"
    script:
        os.path.join(workflow.basedir, "scripts/time_parser/time_compiler.py")