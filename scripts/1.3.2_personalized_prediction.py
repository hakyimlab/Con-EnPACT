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

from pyliftover import LiftOver
import h5py

import pickle as pkl
import sqlite3

import glob

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
intermediates_dir_predictDB_standard = os.path.join(project_dir, "intermediates",
                        "predictDB_standard")
intermediates_dir = os.path.join(project_dir, "intermediates",
                        "personalized_prediction")
os.makedirs(intermediates_dir, exist_ok=True)

context = parameters["general_parameters"]["context"]

training_parameters = parameters["train_enpact"]

run_sbatch = sys.argv[2]


#################################################################
# 1.) Set up individuals and features to predict on
#################################################################

print("Setting up individuals and features")

count_table_file = parameters["generate_enpact_training_data"]["input_files"]["normalized_expression_data"]

count_table = pd.read_csv(count_table_file, sep="\t", index_col=0)

individuals = list(count_table.columns)
features = list(count_table.index)

# Features determined from the stadard predictDB well-predicted set

path_to_filtered_db = os.path.join(intermediates_dir_predictDB_standard, "predictDB",
                                   "filtered_db", "predict_db_Model_training_filtered.db")

def read_from_sqlite(db, query):
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute(query)
    return c.fetchall()

query = "SELECT gene FROM extra"

features = read_from_sqlite(path_to_filtered_db, query)
features_list = [x[0] for x in features]

if os.path.exists(os.path.join(intermediates_dir, "features.txt")):
    print("Sampled features already exist, not overwriting")

    max_features = parameters["personalized_predictions"]["max_features"]
    if len(features_list) > max_features:
        features_list_sampled = random.sample(features_list, max_features)

    with open(os.path.join(intermediates_dir, "features.txt"), "w") as f:
        for feature in features_list_sampled:
            f.write(feature + "\n")

# We also want to reformat the selected features for Epigenome prediction

if parameters["personalized_predictions"]["liftover"]:

    ori_features = []
    lo_features = []

    if parameters["personalized_predictions"]["liftover_target"] == "hg38":
        liftover = LiftOver('hg19', 'hg38')
    elif parameters["personalized_predictions"]["liftover_target"] == "hg19":
        liftover = LiftOver('hg38', 'hg19')

with open(os.path.join(intermediates_dir, "interval_list.txt"), "w") as f:
    for feature in features_list_sampled:

        chrom = feature.split("_")[0]

        feature_split = feature.split("_")

        if parameters["personalized_predictions"]["liftover"]:

            lifted_over_feature_start = liftover.convert_coordinate(feature_split[0], int(feature_split[1]))
            lifted_over_feature_end = liftover.convert_coordinate(feature_split[0], int(feature_split[2]))

            if len(lifted_over_feature_start) == 0 or len(lifted_over_feature_end) == 0:
                print("Liftover failed for feature: ", feature)
                continue

            lifted_over_feature_list = [chrom, lifted_over_feature_start[0][1], lifted_over_feature_end[0][1]]

            center_coord = (lifted_over_feature_list[1] + lifted_over_feature_list[2]) // 2

            ori_features.append(feature)
            lo_features.append("_".join(lifted_over_feature_list))

        else:
            center_coord = (int(feature.split("_")[1]) + int(feature.split("_")[2])) // 2

        f.write(chrom + "_" + str(center_coord) +"_" + str(center_coord+2)+"\n")

if parameters["personalized_predictions"]["liftover"]:
    lo_mapping_df = pd.DataFrame({"ori_features": ori_features, "lo_features": lo_features})
    lo_mapping_df.to_csv(os.path.join(intermediates_dir, "lo_mapping.csv"), index=False)

with open(os.path.join(intermediates_dir, "individuals.txt"), "w") as f:
    for individual in individuals:
        f.write(individual.split("_")[0] + "\n")


#################################################################
# 2.) Set up personalized epigenome predictions
#################################################################
        
print("Setting up personalized epigenome predictions")


path_to_epigenome_prediction_pipeline = parameters["personalized_predictions"]["path_to_epigenome_prediction_pipeline"]
path_to_epigenome_config = parameters["personalized_predictions"]["path_to_epigenome_config"]

loaded_epigenome_config = json.load(open(path_to_epigenome_config))

loaded_epigenome_config["interval_list_file"] = os.path.join(intermediates_dir, "interval_list.txt")
loaded_epigenome_config["individuals"] = os.path.join(intermediates_dir, "individuals.txt")
loaded_epigenome_config["n_individuals"] = -1
loaded_epigenome_config["project_dir"] = intermediates_dir
loaded_epigenome_config["prediction_data_name"] = "predicted_epigenome"
loaded_epigenome_config["prediction_id"] = context

os.makedirs(os.path.join(intermediates_dir,'metadata'), exist_ok=True)
loaded_epigenome_config["metadata_dir"] = os.path.join(intermediates_dir,'metadata')

epigenome_config_parameters = parameters["personalized_predictions"]["epigenome_config_parameters"]
for key in epigenome_config_parameters:
    loaded_epigenome_config[key] = epigenome_config_parameters[key]

path_to_vcf = parameters["personalized_predictions"]["path_to_vcf"]

vcf_files = glob.glob(path_to_vcf + "/*.vcf.gz")

for chrom in loaded_epigenome_config["vcf_files"]["files"].keys():
    for vcf_file in vcf_files:
        if chrom in vcf_file:
            loaded_epigenome_config["vcf_files"]["files"][chrom] = vcf_file


with open(os.path.join(intermediates_dir, "epigenome_config.json"), "w") as f:
    json.dump(loaded_epigenome_config, f, indent=4)



#################################################################
# 2.) Collect personalized Epigenome predictions into txt files for EnPACT
#################################################################

print("Collecting personalized Epigenome predictions")

#################################################################
# 3.) Analyze personalized prediction accuracy of EnPACT vs PredictDB
#################################################################

print("Analyzing personalized prediction accuracy")