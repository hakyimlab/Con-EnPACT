import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd

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
                        "evaluate_training")
os.makedirs(intermediates_dir, exist_ok=True)

context = parameters["general_parameters"]["context"]

training_parameters = parameters["train_enpact"]


#################################################################
# 1.) Run evaluation script
#################################################################

print("Running training script")

evaluation_parameters = training_parameters["evaluation_parameters"]

script_directory = os.path.dirname(os.path.abspath(sys.argv[0])) 
print(script_directory)

if training_parameters["model_type"] == "elastic_net":

    training_script = os.path.join(script_directory,"evaluate_training_elastic_net.R")
    train_data_dir = intermediates_dir_generate_data
    rds_file = os.path.join(intermediates_dir_train_enpact,
                            f"trained_enpact_eln_{context}.linear.rds")
    gene_annotations = parameters["generate_enpact_training_data"]["input_files"]["gene_annotations"]
    epigenome_reference_track = evaluation_parameters["epigenome_reference_track"]
    output_dir = intermediates_dir
    color_palette = parameters["general_parameters"]["color_palette"]

    run_command = [
        "Rscript",
        training_script,
        f"--train_data_dir={train_data_dir}",
        f"--rds_file={rds_file}",
        f"--gene_annotations={gene_annotations}",
        f"--context={context}",
        f"--epigenome_reference_track={epigenome_reference_track}",
        f"--output_dir={output_dir}",
        f"--color_palette={color_palette}"

    ]

    print(" ".join(run_command))

    subprocess.run(run_command)