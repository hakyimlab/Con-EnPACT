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
intermediates_dir = os.path.join(project_dir, "intermediates",
                       "train_enpact")
os.makedirs(intermediates_dir, exist_ok=True)

context = parameters["general_parameters"]["context"]

training_parameters = parameters["train_enpact"]

script_directory = os.path.dirname(os.path.abspath(sys.argv[0])) 
print(script_directory)

#################################################################
# 1.) Run training script
#################################################################

print("Running training script")

if training_parameters["model_type"] == "elastic_net":
    training_script = os.path.join(script_directory,"train_EnPACT_elastic_net.R")
    
    train_data_X = os.path.join(intermediates_dir_generate_data,
                                f"epigenome_train.txt")
    train_data_y = os.path.join(intermediates_dir_generate_data,
                                f"train_{context}_mean_expression.tsv")
    rds_file = os.path.join(intermediates_dir,
                            f"trained_enpact_eln_{context}")

    run_command = [
        "Rscript",
        training_script,
        f"--train_data_file_X={train_data_X}",
        f"--train_data_file_y={train_data_y}",
        f"--rds_file={rds_file}"
    ]

    print(" ".join(run_command))

    subprocess.run(run_command)