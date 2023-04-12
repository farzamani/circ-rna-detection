import os

# configfile: "config.yaml"
SAMPLES = ["A_250k", "A_500k", "A_1m", "A_2m"]
SAMPLES_TEST = ["A_500k"]


rule all:
    input:
        #expand("results/{sample}/CircDBG/time.txt", sample=SAMPLES),
        #expand("results/{sample}/CIRCexplorer2/time.txt", sample=SAMPLES),
        #expand("results/{sample}/CircMarker/time.txt", sample=SAMPLES),
        #expand("results/{sample}/CirComPara2/time.txt", sample=SAMPLES_TEST),
        #expand("results/{sample}/circRNA_finder/time.txt", sample=SAMPLES),
        #expand("results/{sample}/circRNAFinder/time.txt", sample=SAMPLES),
        #expand("results/{sample}/CIRI2/time.txt", sample=SAMPLES),
        #expand("results/{sample}/DCC/time.txt", sample=SAMPLES),
        #expand("results/{sample}/find_circ/time.txt", sample=SAMPLES),
        #expand("results/{sample}/MapSplice/time.txt", sample=SAMPLES),
        #expand("results/{sample}/segemehl/time.txt", sample=SAMPLES),
        expand("results/{sample}/time.csv", sample=SAMPLES)
        

rule trim:
    input:
        reads1 = "reads/{sample}_1.fq.gz",
        reads2 = "reads/{sample}_2.fq.gz"
    output:
        trimmed1 = "trimmed/{sample}_1_val_1.fq.gz",
        trimmed2 = "trimmed/{sample}_2_val_2.fq.gz",
        report1 = temp("trimmed/{sample}_1.fq.gz_trimming_report.txt"),
        report2 = temp("trimmed/{sample}_2.fq.gz_trimming_report.txt")
    conda:
        os.path.join(workflow.basedir, "envs/snakemake.yaml")
    resources:
        mem_mb = 4000,
        time = 720
    shell:
        "trim_galore --paired {input.reads1} {input.reads2} -o trimmed"


include: "rules/CIRI2.smk"
include: "rules/CircMarker.smk"
include: "rules/CircDBG.smk"
include: "rules/segemehl.smk"
include: "rules/DCC.smk"
include: "rules/CIRCexplorer2.smk"
include: "rules/find_circ.smk"
include: "rules/MapSplice.smk"
include: "rules/circRNAFinder.smk"
include: "rules/circRNA_finder.smk"
include: "rules/CirComPara2.smk"
include: "rules/parse_time.smk"