import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd

import json

import kipoiseq

import multiprocessing as mp

import pyBigWig 
import os,sys

import h5py

import pickle as pkl

import epigenome_utils
import logging_utils

#################################################################
# 0.) Load required input variables, and create directories
#################################################################

print("Loading input variables")

json_file = sys.argv[1]
parameters = json.load(open(json_file))

project_dir = parameters["general_parameters"]["project_directory"]

intermediates_dir = os.path.join(project_dir, "intermediates",
                       "generate_enpact_training_data")
os.makedirs(intermediates_dir, exist_ok=True)

context = parameters["general_parameters"]["context"]

input_files = parameters["generate_enpact_training_data"]["input_files"]

chrom_split = parameters["generate_enpact_training_data"]["chromosomes_for_training"]

epigenome_parameters = parameters["generate_enpact_training_data"]["reference_epigenome"]


#################################################################
# 1.) Establish common genes between expression and annotations
#################################################################

print("Establishing common genes between expression and annotations")

gene_annotations = pd.read_csv(input_files["gene_annotations"], index_col=0)
normalized_expression_data = pd.read_csv(input_files["normalized_expression_data"], sep="\t", index_col=0)

common_genes = list(set(gene_annotations.index).intersection(set(normalized_expression_data.index)))
with open(os.path.join(project_dir, "intermediates",
                       "generate_enpact_training_data",
                       "common_genes.txt"), "w") as f:
    for gene in common_genes:
        f.write(gene + "\n")


#################################################################
# 2.) Create train, test, validation splits
#################################################################

print("Creating train, test, validation splits")

gene_mapping = {
    "valid":set(chrom_split["validation"].split(",")),
    "test":set(chrom_split["test"].split(",")),
    "train":set(chrom_split["train"].split(","))
}

partitioned_genes = {
    "train": [],
    "valid": [],
    "test": []
}

training_gene_metadata_file = os.path.join(intermediates_dir,"training_gene_metadata.tsv")

with open(training_gene_metadata_file, "w") as region_file:
    for gene in common_genes:

        gene_chr = "chr"+gene_annotations[gene_annotations.index == gene]["chromosome_name"].values[0]
        if gene_chr in ["chrX","chrY","chrM"]:
            continue
        if gene_chr in gene_mapping["valid"]:
            cur_group = "valid"
            partitioned_genes["valid"].append(gene)
        elif gene_chr in gene_mapping["test"]:
            cur_group = "test"
            partitioned_genes["test"].append(gene)
        elif gene_chr in gene_mapping["train"]:
            cur_group = "train"
            partitioned_genes["train"].append(gene)
        else:
            print(f"Gene {gene_chr} not found in any partition")

        tss_site = gene_annotations[gene_annotations.index == gene]["transcription_start_site"].values[0]
 
        out = [
            gene_chr,
            str(tss_site),
            str(tss_site+2),
            cur_group,
            gene
        ]

        region_file.write("\t".join(out)+"\n")


#################################################################
# 3.) Split expression data into training splits
#################################################################
        
print("Splitting expression data into training splits")

for dset in ["train", "valid", "test"]:
    cur_table = normalized_expression_data[normalized_expression_data.index.isin(partitioned_genes[dset])]
    print(dset, cur_table.shape)
    cur_table.sort_index(inplace=True)

    o_file = os.path.join(intermediates_dir, f"{dset}_{context}_mean_expression.tsv")
    cur_table.mean(axis=1).to_csv(o_file, sep="\t", header=False)


#################################################################
# 4.) Query reference epigenome for training splits
#################################################################

print("Querying reference epigenome for training splits")

ref_epi_dir = epigenome_parameters["reference_epigenome_path"]
num_bins = epigenome_parameters["num_bins"]
tracks = epigenome_parameters["tracks"]

for dset in ["train", "valid", "test"]:
    if os.path.exists(os.path.join(intermediates_dir, f"epigenome_{dset}.txt")):
        os.remove(os.path.join(intermediates_dir, f"epigenome_{dset}.txt"))

genes_set = set()

with open(training_gene_metadata_file, "r") as rf:
    for line in rf:
        p_line = line.strip().split("\t")

        gene = p_line[-1]

        if gene in genes_set:
            print(gene)
            sys.exit()
        else:
            genes_set.add(gene)

        # print(gene)

        if p_line[0] in gene_mapping["valid"]:
            cur_group = "valid"
        elif p_line[0] in gene_mapping["test"]:
            cur_group = "test"
        elif p_line[0] in gene_mapping["train"]:
            cur_group = "train"

        cur_query = epigenome_utils.query_epigenome(p_line[0].lstrip("chr"),
                                     int(p_line[1])+1, 
                                     ref_epi_dir,
                                     num_bins=num_bins, 
                                     tracks=tracks)
        query_vec = [gene]+[str(x) for x  in list(cur_query.mean(axis=0))]

        with open(os.path.join(intermediates_dir, f"epigenome_{cur_group}.txt"), "a") as f:
            f.write("\t".join(query_vec)+"\n")
