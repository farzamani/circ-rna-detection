rule detect_segemehl:
    input:
        reads1 = rules.trim.output.trimmed1,
        reads2 = rules.trim.output.trimmed2,
    output:
        sam = "results/{sample}/segemehl/{sample}.sam",
        time = "results/{sample}/segemehl/time.txt"
    params:
        fa = "ref/reference/hg19.fa",
        index = "ref/reference/hg19.idx",
        time = "results/{sample}/segemehl/time.txt",
        group = "{sample}",
    threads: 
        1
    resources:
        mem_mb = 42000,
        time = 720
    shell:
        "/usr/bin/time -v -o {params.time} bash -c 'segemehl.x -i {params.index} -d {params.fa} -t {threads} -g {params.group} -q {input.reads1} -p {input.reads2} > {output.sam}'"