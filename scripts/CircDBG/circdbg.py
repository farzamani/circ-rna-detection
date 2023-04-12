import os

def generate_config(reads1, reads2, output_file, fa, gtf):
    with open(output_file, "w") as f:
        f.write("[General]\n")
        f.write(f"Reference=../../../{fa}\n")
        f.write(f"GTF=../../../{gtf}\n")
        f.write(f"Reads1=../../../{reads1}\n")
        f.write(f"Reads2=../../../{reads2}\n")
        f.write("[StdCircularRNAInfo]\n")
        f.write("CircRNADb=\n")
        f.write("CircRNAdbResult=\n")
        f.write("Tissue=\n")
        f.write("[Parameter]\n")
        f.write("KmerLen=15\n")
        f.write("ReadLen=101\n")
        f.write("KmerRatio=60\n")
        f.write("ThreadsNum=4\n")
        f.write("DoDetection=1\n")
        f.write("DoComparison=0\n")
        f.write("MaxSupportNum=999\n")
        f.write("[Comparison]\n")
        f.write("MyResult=\n")
        f.write("CiriResult=\n")
        f.write("CIRCExplorerResult=\n")
        f.write("FindCircResult=\n")
        f.write("CircRNAFinderResult=\n")
        f.write("CircMarkerResult=\n")
        f.write("SimulatorBenchmark=")

    return output_file


print(snakemake.input.reads1)
print(snakemake.input["reads1"])
print(snakemake.input["reads2"])
print(snakemake.output["config"])
print(snakemake.params["fa"])
print(snakemake.params["gtf"])

generate_config(snakemake.input["reads1"], snakemake.input["reads2"], snakemake.output["config"], snakemake.params["fa"], snakemake.params["gtf"])
