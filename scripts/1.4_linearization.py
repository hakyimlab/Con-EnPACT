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
os.makedirs(os.path.join(intermediates_dir,"GEUVADIS_predictions_mean"), exist_ok=True)

context = parameters["general_parameters"]["context"]

linearization_parameters = parameters["linearization_parameters"]


#################################################################
# 1.) Prepare personalized epigenome for EnPACT inference
#################################################################

if linearization_parameters["linearization_dataset"] == "GEUVADIS":

    with open(linearization_parameters["individuals"], "r") as inds_f:
        inds = inds_f.read().split()

enformer_geuvadis_pred_folder = linearization_parameters["epigenome_pred_dir"]

with h5py.File(os.path.join(enformer_geuvadis_pred_folder, inds[0]+".h5"), "r") as f:
    genes_dsets = list(f.keys())
    genes_dsets.sort()

for ind in inds:
    if os.path.exists(os.path.join(intermediates_dir,"GEUVADIS_predictions_mean",ind+".txt")):
        continue
    expression_array = np.zeros((len(genes_dsets),5313))
    if os.path.exists(os.path.join(enformer_geuvadis_pred_folder, ind+".h5")):
        with h5py.File(os.path.join(enformer_geuvadis_pred_folder, ind+".h5"), "r") as f:
            for i, gene in enumerate(genes_dsets):
                # print(i)
                # print(np.mean(f[gene][:,:],axis=0).shape)
                expression_array[i,:] = np.mean(f[gene][:,:],axis=0)
            np.savetxt(os.path.join("/beagle3/haky/users/saideep/projects/aracena_modeling/linearization",ind+".txt"),expression_array)
    else:
        print(ind)


#################################################################
# 2.)  Predict for linearization
#################################################################
