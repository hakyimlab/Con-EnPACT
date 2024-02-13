import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import h5py

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

linearization_parameters = parameters["linearization_parameters"]


#################################################################
# 1.) Prepare personalized epigenome for EnPACT inference
#################################################################

if linearization_parameters["linearization_dataset"] == "GEUVADIS":

    with open("/beagle3/haky/users/saideep/github_repos/Daily-Blog-Sai/posts/2023-11-16-linearization/individuals.txt", "r") as inds_f:
        inds = inds_f.read().split()
