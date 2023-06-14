import os
import csv
import datetime
import pandas as pd

replicates = snakemake.params.replicate

def compile_files(infile):
    df = pd.read_csv(infile)

    for replicate in replicates:
        df["replicate"] = replicate




for path in snakemake.input.time:
    append_to_csv(snakemake.output.time, parse_file(path))