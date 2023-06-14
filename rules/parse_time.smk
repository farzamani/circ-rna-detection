ALL_TOOLS = [
    "CircDBG", "CIRCexplorer2", "CircMarker", "CirComPara2",
    "circRNA_finder", "circRNAFinder", "CircSplice", "CIRI2", 
    "DCC", "find_circ", "segemehl"
]

TOOLS = [
    "CircDBG", 
    "CircMarker", 
    "circRNAFinder",
    "MapSplice"
]

# REPLICATE = [0, 1, 2, 3, 4]
REPLICATE = [1]

rule parse_time:
    input:
        time = expand("results/{{sample}}.{{replicate}}/{tools}/time.txt", tools=ALL_TOOLS)
    output:
        time = "results/{sample}.{replicate}/time.csv"
    script:
        os.path.join(workflow.basedir, "scripts/time_parser/time_parser.py")


rule compile_time:
    input:
        time = expand("results/{{sample}}.{replicate}/time.csv", replicate=REPLICATE)
    output:
        time = "results/{sample}/time.csv"
    params:
        sample = "{sample}",
        tools = ALL_TOOLS,
        replicate = REPLICATE
    script:
        os.path.join(workflow.basedir, "scripts/time_parser/replicate_compiler.py")