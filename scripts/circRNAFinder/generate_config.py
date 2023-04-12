import os

def generate_config_find(sample_name, reads, output_file, bowtie, transcriptome, fa):
    with open(output_file, "w") as f:
        f.write(f"[lib_{sample_name}]\n")
        f.write(f"sample_name = {sample_name}\n")
        f.write(f"reads_file = {reads}\n")
        f.write(f"genome_index = {bowtie}\n")
        f.write(f"trans_index = {transcriptome}\n")
        f.write(f"genome_seq = {fa}\n")
    return output_file

def generate_config_anno(sample_name, output_file, bowtie, fa, gtf_build):
    with open(output_file, "w") as f:
        f.write(f"[group_{sample_name}]\n")
        f.write(f"sample_names = {sample_name}\n")
        f.write(f"genome_index = {bowtie}\n")
        f.write(f"genome_seq = {fa}\n")
        f.write(f"gene_bed = {gtf_build}\n")
        f.write(f"id_prefix = {sample_name}\n")
    return output_file

print("Generating config find")
generate_config_find(snakemake.params.sample_name, 
                        snakemake.input.reads, 
                        snakemake.output.config_find, 
                        snakemake.params.bowtie, 
                        snakemake.params.transcriptome, 
                        snakemake.params.fa)

print("Generating config anno")
generate_config_anno(snakemake.params.sample_name, 
                        snakemake.output.config_anno, 
                        snakemake.params.bowtie, 
                        snakemake.params.fa, 
                        snakemake.params.gtf_build)