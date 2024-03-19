import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import h5py
import time

import subprocess

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

intermediates_dir_generate_data = os.path.join(project_dir, "intermediates",
                       "generate_enpact_training_data")
intermediates_dir_train_enpact = os.path.join(project_dir, "intermediates",
                       "train_enpact")
intermediates_dir = os.path.join(project_dir, "intermediates",
                        "linearization")
os.makedirs(intermediates_dir, exist_ok=True)

context = parameters["general_parameters"]["context"]

linearization_parameters = parameters["linearization"]

script_directory = os.path.dirname(os.path.abspath(sys.argv[0])) 
print(script_directory)

#################################################################
# 1.) Prepare personalized epigenome for EnPACT predictions
#################################################################

print("Preparing personalized epigenome for EnPACT predictions")

def create_inputs_for_ind(ind, intermediates_dir, epigenome_pred_folder, genes_dsets, num_bins, num_tracks):
    print(mp.current_process(), " individual: ", ind)
    start_time = time.time()
    expression_array = np.zeros((len(genes_dsets),num_tracks))
    if os.path.exists(os.path.join(intermediates_dir,"personalized_epigenome_predictions_mean_"+ind+".txt")):
        print(ind," already exists")
        return
    if os.path.exists(os.path.join(epigenome_pred_folder, ind+".h5")):
        with h5py.File(os.path.join(epigenome_pred_folder, ind+".h5"), "r") as f:
            for i, gene in enumerate(genes_dsets):
                expression_array_bins = f[gene].shape[0]
                central_bin_start_ind = find_central_bin_start_ind(num_bins, expression_array_bins)
                expression_array[i,:] = np.mean(f[gene][central_bin_start_ind:central_bin_start_ind+num_bins,:],axis=0)
            np.savetxt(os.path.join(intermediates_dir,"personalized_epigenome_predictions_mean_"+ind+".txt"),expression_array)
    
    end_time = time.time()
    print("Time taken for ", ind, " is ", end_time-start_time)

def find_central_bin_start_ind(mid_len, all_bins):
    return ((all_bins-mid_len)//2)


with open(linearization_parameters["individuals"], "r") as inds_f:
    inds = inds_f.read().split()

epigenome_pred_folder = linearization_parameters["epigenome_pred_dir"]

num_tracks = parameters["generate_enpact_training_data"]["reference_epigenome"]["num_tracks"]
num_bins = parameters["generate_enpact_training_data"]["reference_epigenome"]["num_bins"]

with h5py.File(os.path.join(epigenome_pred_folder, inds[0]+".h5"), "r") as f:
    genes_dsets = list(f.keys())
    genes_dsets.sort()

args = []
for ind in inds:
    print(ind)
    args.append((ind, intermediates_dir, epigenome_pred_folder, genes_dsets, num_bins, num_tracks))
cpu_count = mp.cpu_count()
with mp.Pool(cpu_count) as pool:
    print(cpu_count)
    pool.starmap(create_inputs_for_ind, args)


#################################################################
# 2.)  Make personalized EnPACT predictions for linearization
#################################################################


def make_personalized_enpact_predictions(ind, intermediates_dir, model_path):
    if os.path.exists(os.path.join(intermediates_dir,f"enpact_personalized_predictions_{ind}.txt")):
        print(ind," already exists")
        return

    com = [
        "Rscript",
        "--vanilla",
        os.path.join(script_directory,"predict_for_linearization_elastic_net.R"),
        f"--individual={ind}",
        f"--intermediates_dir={intermediates_dir}",
        f"--trained_model={model_path}"
    ]

    print(" ".join(com))

    subprocess.run(com)

model_path = os.path.join(intermediates_dir_train_enpact,
                        f"trained_enpact_eln_{context}.linear.rds")

args = []
for ind in inds: 
    args.append((ind, intermediates_dir, model_path))
with mp.Pool(processes=mp.cpu_count()) as pool:
    pool.starmap(make_personalized_enpact_predictions, args)


#################################################################
# 3.)  Collect personalized EnPACT predictions for linearization
#################################################################

mode = "Flu"

with open("/beagle3/haky/users/saideep/github_repos/Daily-Blog-Sai/posts/2023-11-16-linearization/individuals.txt", "r") as inds_f:
    inds = inds_f.read().split()


paste_com = [
    "paste",
    "-d'\t'"
]

inds_include = []
for ind in inds:
    if os.path.exists(os.path.join(intermediates_dir,f"enpact_personalized_predictions_{ind}.txt")):
        paste_com.append(os.path.join(intermediates_dir,f"enpact_personalized_predictions_{ind}.txt"))
        inds_include.append(ind)

paste_com.append("|")
paste_com.append("tail")
paste_com.append("-n")
paste_com.append("+2")

paste_com.append(">")
paste_com.append(os.path.join(intermediates_dir, "enpact_personalized_predictions_combined.txt"))

print(" ".join(paste_com))

os.system(" ".join(paste_com))


import pandas as pd

combined_expression_table = pd.read_csv(os.path.join(intermediates_dir, "enpact_personalized_predictions_combined.txt"), delimiter="\t", header=None)

combined_expression_table.columns = inds_include

with h5py.File(os.path.join(epigenome_pred_folder, inds_include[0]+".h5"), "r") as f:
    genes_dsets = list(f.keys())
    genes_dsets.sort()

# Create mapping of gene regions to ensemble gene ids
gene_region_mapping = {}
c=0
with open(parameters["generate_enpact_training_data"]["input_files"]["gene_annotations"], "r") as genes_f:
    for line in genes_f:
        if c==0:
            c+=1
            continue
        tss_region = line.strip().split(",")[9].split("_")
        tss_plus_one = "_".join([tss_region[0],str(int(tss_region[1])-1),str(int(tss_region[2])-1)])
        gene_region_mapping[tss_plus_one] = line.split(",")[0]


# Filter gene list for genes not in the gene annotation file
include_gene = []
gene_ids = []
for gene in genes_dsets:
    cur_gene_region = gene.rstrip("_predictions")
    if cur_gene_region not in gene_region_mapping:
        print(cur_gene_region)
        include_gene.append(False)
    else:
        gene_ids.append(gene_region_mapping[cur_gene_region])
        include_gene.append(True)

combined_expression_table_genes_in_anno = combined_expression_table[include_gene]
combined_expression_table_genes_in_anno.index = gene_ids

combined_expression_table_genes_in_anno = combined_expression_table_genes_in_anno.dropna(axis=1)

combined_expression_table_genes_in_anno.to_csv(os.path.join(intermediates_dir, "enpact_personalized_predictions_combined_labeled_filtered.txt"),
                                                index=True, index_label="NAME", sep="\t")


#################################################################
# 4.)  Build PredictDB sbatch script to create linearized EnPACT expression models
#################################################################

with open(os.path.join(script_directory,"run_predictDB.sbatch"), "r") as f:
    predictDB_script = f.read()

    predictDB_script = predictDB_script.replace("JOBNAME", "EnPACT_Linearization")

    predictDB_script = predictDB_script.replace("PREDICTDB_REPO", linearization_parameters["path_to_predictDB_pipeline"])
    predictDB_script = predictDB_script.replace("EXPRESSION", os.path.join(intermediates_dir, "enpact_personalized_predictions_combined_labeled_filtered.txt"))
    predictDB_script = predictDB_script.replace("OUTPUT_DIR", os.path.join(intermediates_dir,"predictDB"))
    predictDB_script = predictDB_script.replace("GENE_ANNOTATION", linearization_parameters["gene_annotation"])
    predictDB_script = predictDB_script.replace("SNP_ANNOTATION", linearization_parameters["snp_annotation"])
    predictDB_script = predictDB_script.replace("GENOTYPE", linearization_parameters["genotype_data"])


    os.makedirs(os.path.join(intermediates_dir, "predictDB"), exist_ok=True)
    with open(os.path.join(intermediates_dir, "predictDB","run_predictDB.sbatch"), "w") as rp:
        rp.write(predictDB_script)



    subprocess.run(["sbatch", os.path.join(intermediates_dir, "predictDB","run_predictDB.sbatch")],
                   cwd= os.path.join(intermediates_dir, "predictDB"))