import os

def generate_meta(sample_name, reads1, reads2, output_file):
    with open(output_file, "w") as f:
        f.write("file,sample\n")
        f.write(f"../../../{reads1},{sample_name}\n")
        f.write(f"../../../{reads2},{sample_name}\n")
    return output_file


def generate_vars(output_file, meta, fa, gtf, segemehl, bwa, bowtie2, bowtie, star, hisat2, genepred, threads):
    with open(output_file, "w") as f:
        f.write(f"META = '{meta}'\n")
        f.write(f"GENOME_FASTA = '../../../{fa}'\n")
        f.write(f"ANNOTATION = '../../../{gtf}'\n")
        f.write(f"SEGEMEHL_INDEX = '../../../{segemehl}'\n")
        f.write(f"BWA_INDEX = '../../../{bwa}'\n")
        f.write(f"BOWTIE2_INDEX = '../../../{bowtie2}'\n")
        f.write(f"BOWTIE_INDEX = '../../../{bowtie}'\n")
        f.write(f"STAR_INDEX = '../../../{star}'\n")
        f.write(f"GENOME_INDEX = '../../../{hisat2}'\n")
        f.write(f"GENEPRED = '../../../{genepred}'\n")
        f.write("BYPASS = 'linear'\n")
        f.write(f"CPUS = '../../../{threads}'\n")
    return output_file


print("Generating meta")
generate_meta(snakemake.params.sample_name, 
                snakemake.input.reads1,
                snakemake.input.reads2, 
                snakemake.output.meta_csv)

generate_vars(snakemake.output.vars_py,
                snakemake.params.meta,
                snakemake.params.fa,
                snakemake.params.gtf,
                snakemake.params.segemehl,
                snakemake.params.bwa,
                snakemake.params.bowtie2,
                snakemake.params.bowtie,
                snakemake.params.star,
                snakemake.params.hisat2,
                snakemake.params.genepred,
                snakemake.threads)