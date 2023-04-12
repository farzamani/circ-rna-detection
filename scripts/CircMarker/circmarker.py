import os

def generate_config(reads1, reads2, output_file, fa, gtf, bwa):
    with open(output_file, "w") as f:
        f.write("[General]\n")
        f.write(f"Reference=../../../{fa}\n")
        f.write(f"GTF=../../../{gtf}\n")
        f.write(f"Reads1=../../../{reads1}\n")
        f.write(f"Reads2=../../../{reads2}\n")
        f.write("[Parameters]\n")
        f.write("MinSupportReads=1\n")
        f.write("MaxSupportReads=999\n")
        f.write("ReadsLen=101\n")
        f.write("KmerRatio=30\n")
        f.write("[Results]\n")
        f.write(f"CurResult=\n")
        f.write(f"SimulationResult=\n")
        f.write("CIRIResult=\n")
        f.write("CIRCexplorerResult=\n")
        f.write("FindCircResult=\n")
        f.write("CircBaseResult=\n")
        f.write("[HitCompare]\n")
        f.write("CIRIHit=\n")
        f.write("CurHit=\n")
        f.write("FindCircHit=\n")
        f.write("CircExplorerHit=\n")
        f.write("CurChromIndex=1\n")
        f.write("[Mapping]\n")
        f.write(f"BWA=../../../{bwa}\n")

    return output_file


print(snakemake.input.reads1)
print(snakemake.input["reads1"])
print(snakemake.input["reads2"])
print(snakemake.output["config"])
print(snakemake.params["fa"])
print(snakemake.params["gtf"])

generate_config(snakemake.input["reads1"], snakemake.input["reads2"], snakemake.output["config"], snakemake.params["fa"], snakemake.params["gtf"], snakemake.params.bwa)
