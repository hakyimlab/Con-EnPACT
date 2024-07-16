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

import enpact_utils
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
# 1.) Check predicted epigenome for each linearization dataset
#################################################################

print("Preparing features to linearize EnPACT model on")

linearization_datasets = linearization_parameters["linearization_datasets"]

# for ld in linearization_datasets.keys():


#TODO Add in checks here
#TODO Need to check that all the files exist, including mean summarized enformer predictions

#################################################################
# 2.) Make EnPACT predictions from personalized epigenome
#################################################################

print("Making EnPACT predictions")

window_size = parameters["generate_enpact_training_data"]["reference_epigenome"]["num_bins"]
context = parameters["general_parameters"]["context"]

model_path = os.path.join(intermediates_dir_train_enpact,
                        f"trained_enpact_eln_{context}.linear.rds")

for ld in linearization_datasets.keys():
    epigenome_pred_dir = linearization_datasets[ld]["epigenome_pred_dir"]
    inds_file = linearization_datasets[ld]["individuals"]
    with open(inds_file, "r") as inds_f:
        inds = inds_f.read().split()

    cur_lin_dir = os.path.join(intermediates_dir, ld)
    os.makedirs(os.path.join(cur_lin_dir, "enpact_preds"), exist_ok=True)

    pool = mp.Pool(processes=mp.cpu_count())

    arguments_for_starmap = []
    output_enpact_predictions = []
    for ind in inds:
        print(f"{ld}_{ind}")
        cur_enformer_pred = os.path.join(epigenome_pred_dir, ind+f"_{window_size}.txt")
        if not os.path.exists(cur_enformer_pred):
            print(cur_enformer_pred, " does not exist")
            continue
        cur_enpact_pred = os.path.join(cur_lin_dir, "enpact_preds", f"{ind}_{window_size}_enpact.txt")
        arguments_for_starmap.append((cur_enformer_pred, cur_enpact_pred, model_path, script_directory))
    
    output_enpact_predictions = pool.starmap(enpact_utils.make_enpact_predictions, arguments_for_starmap)
        
# sys.exit()

#################################################################
# 3.) Format EnPACT predictions for PredictDB and prepare inputs
#################################################################

print("Formatting EnPACT predictions for PredictDB")

for ld in linearization_datasets.keys():

    cur_lin_dir = os.path.join(intermediates_dir, ld)

    paste_com = [
        "paste",
        "-d'\t'"
    ]

    inds_file = linearization_datasets[ld]["individuals"]
    with open(inds_file, "r") as inds_f:
        inds = inds_f.read().split()

    inds_include = []
    ind_count = 0
    for ind in inds:

        if os.path.exists(os.path.join(cur_lin_dir, "enpact_preds",f"{ind}_{window_size}_enpact.txt")):
            if ind_count == 0:
                paste_com.append(os.path.join(cur_lin_dir, "enpact_preds",f"{ind}_{window_size}_enpact.txt.rownames"))
                ind_count += 1
            paste_com.append(os.path.join(cur_lin_dir, "enpact_preds",f"{ind}_{window_size}_enpact.txt"))
            inds_include.append(ind)

    paste_com.append("|")
    paste_com.append("tail")
    paste_com.append("-n")
    paste_com.append("+2")

    paste_com.append(">")
    paste_com.append(os.path.join(cur_lin_dir, "enpact_preds", f"combined_predictions_body.txt"))

    # print(" ".join(paste_com))

    os.system(" ".join(paste_com))

    with open(os.path.join(cur_lin_dir, "enpact_preds", f"combined_predictions_header.txt"), "w") as o:
        header = ["NAME"]
        for ind in inds_include:
            header.append(ind)

        o.write("\t".join(header) + "\n")

    os.system(f"cat {os.path.join(cur_lin_dir, 'enpact_preds', 'combined_predictions_header.txt')} {os.path.join(cur_lin_dir, 'enpact_preds', 'combined_predictions_body.txt')} > {os.path.join(cur_lin_dir, 'enpact_preds', 'combined_predictions.txt')}")

    # Clean up table formatting as needed

    combined_preds = pd.read_csv(os.path.join(cur_lin_dir, "enpact_preds", "combined_predictions.txt"), sep="\t")

    combined_preds["NAME"] = [x.split(".")[1] for x in list(combined_preds["NAME"])]

    combined_preds.to_csv(os.path.join(cur_lin_dir, "enpact_preds", "combined_predictions.txt"), sep="\t", index=False)

#################################################################
# 4.) Create PredictDB sbatch scripts
#################################################################

print("Creating PredictDB sbatch scripts")

predictDB_template = os.path.join(script_directory, "run_predictDB.sbatch")

for ld in linearization_datasets.keys():

    cur_lin_dir = os.path.join(intermediates_dir, ld)

    predictDB_dir = os.path.join(cur_lin_dir,"predictDB")
    os.makedirs(predictDB_dir, exist_ok=True)

    with open(os.path.join(predictDB_template), "r") as f:
        predictDB_script = f.read()

        predictDB_script = predictDB_script.replace("JOBNAME", "EnPACT_Linearization")

        predictDB_script = predictDB_script.replace("PREDICTDB_REPO", linearization_parameters["path_to_predictDB_pipeline"])
        predictDB_script = predictDB_script.replace("EXPRESSION", os.path.join(cur_lin_dir, "enpact_preds", f"combined_predictions.txt"))
        predictDB_script = predictDB_script.replace("OUTPUT_DIR", predictDB_dir)
        predictDB_script = predictDB_script.replace("GENE_ANNOTATION", linearization_datasets[ld]["gene_annotation_file"])
        predictDB_script = predictDB_script.replace("SNP_ANNOTATION", linearization_datasets[ld]["snp_annotation_file"])
        predictDB_script = predictDB_script.replace("GENOTYPE", linearization_datasets[ld]["genotype_file"])
        predictDB_script = predictDB_script.replace("NFOLDS", "10")


        with open(os.path.join(predictDB_dir,"run_predictDB.sbatch"), "w") as rp:
            rp.write(predictDB_script)

    subprocess.run(["sbatch", os.path.join(predictDB_dir,"run_predictDB.sbatch")], cwd=predictDB_dir)