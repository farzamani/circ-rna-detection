TOOLS = [
    "CircDBG", "CIRCexplorer2", "CircMarker", "CirComPara2",
    "circRNA_finder", "circRNAFinder", "CircSplice", "CIRI2", 
    "DCC", "find_circ", "MapSplice", "segemehl"
]   

rule parse_time:
    input:
        time = expand("results/{{sample}}/{tools}/time.txt", tools=TOOLS)
    output:
        time = "results/{sample}/time.csv"
    script:
        os.path.join(workflow.basedir, "scripts/time_parser/time_parser.py")