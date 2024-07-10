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
import enpact_prediction
import plotting_utils


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

script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))

skip_bool_arg = sys.argv[2]
skip_bool = True
if skip_bool_arg == "False":
    print("skipping")
    skip_bool = False



#################################################################
# 1.) Set up individuals and features to predict on
#################################################################

print("Setting up individuals and features")

count_table_file = parameters["generate_enpact_training_data"]["input_files"]["normalized_expression_data"]

count_table = pd.read_csv(count_table_file, sep="\t", index_col=0)

individuals = list(count_table.columns)
full_features_list = list(count_table.index)

# Features determined from the standard predictDB well-predicted set

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

max_features = parameters["personalized_predictions"]["max_features"]

# If sampled features already exist, skip
if os.path.exists(os.path.join(intermediates_dir, "features.txt")):
    print("Sampled features already exist, not overwriting")
else:
    if len(features_list) > max_features:
        features_list_sampled = random.sample(features_list, max_features)
    else:
        features_list_sampled = features_list

    non_predictdb_features = set(full_features_list) - set(features_list)

    with open(os.path.join(intermediates_dir, "features.txt"), "w") as f:
        for feature in features_list_sampled:
            f.write(feature + "\n")
    
# Non-predictDB features can also be predicted on
if os.path.exists(os.path.join(intermediates_dir, "non_predictdb_features.txt")):
    print("Sampled non_predictdb features already exist, not overwriting")
else:
    non_predictdb_features = set(full_features_list) - set(features_list)
    if len(non_predictdb_features) > max_features:
        non_predictdb_features_sampled = random.sample(non_predictdb_features, max_features)
    else:
        non_predictdb_features_sampled = non_predictdb_features

    with open(os.path.join(intermediates_dir, "non_predictdb_features.txt"), "w") as f:
        for npdb_feature in non_predictdb_features_sampled:
            f.write(npdb_feature + "\n")

# We also want to reformat the selected features for Epigenome prediction

with open(os.path.join(intermediates_dir, "features.txt"), "r") as f:
    features_list_sampled = [x.strip() for x in f]

with open(os.path.join(intermediates_dir, "non_predictdb_features.txt"), "r") as f:
    non_predictdb_features_sampled = [x.strip() for x in f]

if parameters["personalized_predictions"]["liftover"]:

    ori_features = []
    lo_features = []
    interval_coords = []

    if parameters["personalized_predictions"]["liftover_target"] == "hg38":
        liftover = LiftOver('hg19', 'hg38')
    elif parameters["personalized_predictions"]["liftover_target"] == "hg19":
        liftover = LiftOver('hg38', 'hg19')


with open(os.path.join(intermediates_dir, "interval_list.txt"), "w") as f:

    # Loop through the features and write the center coordinates of the features to a file
    feature_lists = [features_list_sampled, non_predictdb_features_sampled]
    for cur_feature_list in feature_lists:
        for feature in cur_feature_list:

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
                lo_features.append("_".join([str(x) for x in lifted_over_feature_list]))
                interval_coords.append("_".join([chrom, str(center_coord), str(center_coord+2)]))


            else:
                center_coord = (int(feature.split("_")[1]) + int(feature.split("_")[2])) // 2

            f.write(chrom + "_" + str(center_coord) +"_" + str(center_coord+2)+"\n")

if parameters["personalized_predictions"]["liftover"]:
    lo_mapping_df = pd.DataFrame({"ori_features": ori_features, "lo_features": lo_features, "interval_coord": interval_coords})
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
loaded_epigenome_config["n_regions"] = -1
loaded_epigenome_config["date"] = parameters["personalized_predictions"]["date"]
loaded_epigenome_config["project_dir"] = intermediates_dir
loaded_epigenome_config["prediction_data_name"] = "predicted_epigenome"
loaded_epigenome_config["prediction_id"] = context
loaded_epigenome_config["batch_individuals"] = -1

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
            skip_bool_vcf = False
            if chrom=="chr1":
                for x in range(10,20):
                    if "chr"+str(x) in vcf_file:
                        skip_bool_vcf = True
                        break
            if chrom=="chr2":
                for x in range(20,30):
                    if "chr"+str(x) in vcf_file:
                        skip_bool_vcf = True
                        break
            if skip_bool_vcf:
                continue

            loaded_epigenome_config["vcf_files"]["files"][chrom] = vcf_file

with open(os.path.join(intermediates_dir, "epigenome_config.json"), "w") as f:
    json.dump(loaded_epigenome_config, f, indent=4)


#################################################################
# 3.) Collect personalized Epigenome predictions into txt files for EnPACT
#################################################################

print("Collecting personalized Epigenome predictions")

date = parameters["personalized_predictions"]["date"]
path_to_predictions = os.path.join(intermediates_dir, "predictions_folder",f"predicted_epigenome_{context}",f"predictions_{date}","enformer_predictions")
collected_preds_dir = os.path.join(intermediates_dir, "collected_predictions")
os.makedirs(collected_preds_dir, exist_ok=True)

window_size = parameters["personalized_predictions"]["num_bins"]
center_bin = 896//2
half_bins = window_size // 2
start_bin = center_bin - half_bins
end_bin = center_bin + half_bins

print(path_to_predictions)

regions = []
with open(os.path.join(intermediates_dir, "interval_list.txt"), "r") as int_list:
    for line in int_list:
        regions.append(line.strip())

# Define function to collect epigenome for given individual
        
def collect_epigenome(ind, regions, path_to_predictions, collected_preds_dir, start_bin, end_bin):

    print("Collecting epigenome for individual: ",ind)

    predictions_dict_haplo1 = {}
    predictions_dict_haplo2 = {}

    for i,region in enumerate(regions):

        pred_file_h1 = os.path.join(path_to_predictions,ind,"haplotype1",region+"_predictions.h5")
        pred_file_h2 = os.path.join(path_to_predictions,ind,"haplotype2",region+"_predictions.h5")

        if not (os.path.exists(pred_file_h1) and os.path.exists(pred_file_h2)):
            continue
        
        try:
            with h5py.File(pred_file_h1,"r") as f:
                pred_h1 = f[region][:]
                pred_h1 = pred_h1[start_bin:end_bin,:]
                predictions_dict_haplo1[region] = list(pred_h1.mean(axis=0))
            with h5py.File(pred_file_h2,"r") as f:
                pred_h2 = f[region][:]
                pred_h2 = pred_h2[start_bin:end_bin,:]
                predictions_dict_haplo2[region] = list(pred_h2.mean(axis=0))
        except:
            print("Error in reading file: ",ind,region)
            continue

    predictions_df = pd.DataFrame(predictions_dict_haplo1).T
    predictions_df.to_csv(os.path.join(collected_preds_dir, ind+"_haplo1.txt"),sep="\t", header=False)

    predictions_df = pd.DataFrame(predictions_dict_haplo2).T
    predictions_df.to_csv(os.path.join(collected_preds_dir, ind+"_haplo2.txt"),sep="\t", header=False)

    print("Finished collecting epigenome for individual: ",ind)

# Compute number of regions and individuals passable

if not skip_bool:
    for ind in individuals:
        skip_count = 0
        for i,region in enumerate(regions):
            pred_file_h1 = os.path.join(path_to_predictions,ind,"haplotype1",region+"_predictions.h5")
            pred_file_h2 = os.path.join(path_to_predictions,ind,"haplotype2",region+"_predictions.h5")

            if not (os.path.exists(pred_file_h1) and os.path.exists(pred_file_h2)):
                skip_count += 1

        print(f"Skipping {skip_count} out of {len(regions)} regions for individual {ind}")

collection_arguments_for_starmap = []

for ind in individuals:
    if skip_bool:
        if os.path.exists(os.path.join(collected_preds_dir, ind+"_haplo1.txt")):
            if os.path.exists(os.path.join(collected_preds_dir, ind+"_haplo2.txt")):
                print("Already exists: ",ind)
                continue
 
    collection_arguments_for_starmap.append((ind, regions, path_to_predictions, collected_preds_dir, start_bin, end_bin))

num_cpus = mp.cpu_count()
print("Using ",num_cpus," CPUs")

with mp.Pool(num_cpus) as pool:
    pool.starmap(collect_epigenome, collection_arguments_for_starmap)


#################################################################
# 4.) Make personalized EnPACT predictions
#################################################################

print("Making personalized EnPACT predictions")

enpact_preds_dir = os.path.join(intermediates_dir, "enpact_predictions")
os.makedirs(enpact_preds_dir, exist_ok=True)

if training_parameters["model_type"] == "elastic_net":
    model_tag = "linear.rds"
elif training_parameters["model_type"] == "logistic":
    model_tag = "logistic.rds"

path_to_trained_model = os.path.join(intermediates_dir_train_enpact, f"trained_enpact_eln_{context}.{model_tag}")

for ind in individuals:
    haplo1_epigenome = os.path.join(collected_preds_dir, ind+"_haplo1.txt")
    haplo1_enpact_pred = os.path.join(enpact_preds_dir, ind+"_haplo1.txt")

    haplo2_epigenome = os.path.join(collected_preds_dir, ind+"_haplo2.txt")
    haplo2_enpact_pred = os.path.join(enpact_preds_dir, ind+"_haplo2.txt")

    if skip_bool:
        if os.path.exists(haplo1_enpact_pred) and os.path.exists(haplo2_enpact_pred):
            print("Already exists: ",ind)
            continue

    enpact_prediction.make_enpact_prediction(haplo1_epigenome, haplo1_enpact_pred, 
                                            path_to_trained_model, 
                                            script_directory)
    
    enpact_prediction.make_enpact_prediction(haplo2_epigenome, haplo2_enpact_pred,
                                            path_to_trained_model,
                                            script_directory)

#################################################################
# 5.) Analyze personalized prediction accuracy of EnPACT vs PredictDB
#################################################################

print("Analyzing personalized prediction accuracy")

pp_analysis_dir = os.path.join(intermediates_dir, "personalized_prediction_accuracy")
os.makedirs(pp_analysis_dir, exist_ok=True)

pp_analysis_plot_dir = os.path.join(pp_analysis_dir, "plots")
os.makedirs(pp_analysis_plot_dir, exist_ok=True)

# Merge the predictions of the two haplotypes within individuals and make a single enpact pred table

if os.path.exists(os.path.join(intermediates_dir, "mean_preds.txt")) and skip_bool:

    mean_preds_df = pd.read_csv(os.path.join(intermediates_dir, "mean_preds.txt"), sep="\t", index_col=0)

else:

    mean_preds_dict = {}

    for ind in individuals:

        haplo1_enpact_pred = os.path.join(enpact_preds_dir, ind+"_haplo1.txt")
        haplo2_enpact_pred = os.path.join(enpact_preds_dir, ind+"_haplo2.txt")

        if not (os.path.exists(haplo1_enpact_pred) and os.path.exists(haplo2_enpact_pred)):
            continue

        haplo1_df = pd.read_csv(haplo1_enpact_pred, sep="\t", header=None, index_col=0)
        haplo2_df = pd.read_csv(haplo2_enpact_pred, sep="\t", header=None, index_col=0)

        # Check if DFs are empty

        if haplo1_df.shape[0] == 0 or haplo2_df.shape[0] == 0:
            continue

        inner_joined_df = haplo1_df.join(haplo2_df, lsuffix="_haplo1", rsuffix="_haplo2")

        inner_joined_df["mean_pred"] = np.mean(inner_joined_df[["1_haplo1", "1_haplo2"]], axis=1)

        if inner_joined_df.shape[0] == 2:
            continue
        mean_preds_dict[ind] = inner_joined_df["mean_pred"]

    mean_preds_df = pd.DataFrame(mean_preds_dict)
    mean_preds_df.to_csv(os.path.join(intermediates_dir, "mean_preds.txt"), sep="\t")

# Filter samples with high proportion of missing values

max_na_proportion = 0.5
mean_preds_df = mean_preds_df.loc[:, mean_preds_df.isna().mean() < max_na_proportion]    

# If lifted over, we need to convert the predictions back to the original coordinates
if parameters["personalized_predictions"]["liftover"]:
    lo_mapping_dict = {}
    with open(os.path.join(intermediates_dir, "lo_mapping.csv"), "r") as f:
        for line in f:
            ori_feature, lo_feature, interval_coord = line.strip().split(",")
            lo_mapping_dict[interval_coord] = ori_feature

for interval in mean_preds_df.index:
    if parameters["personalized_predictions"]["liftover"]:
        mean_preds_df.rename(index={interval: lo_mapping_dict[interval]}, inplace=True)

# Load the predictDB predictions

path_to_predictDB_predictions = os.path.join(intermediates_dir_predictDB_standard,"personalized_predictions.txt")

predictDB_predictions_df = pd.read_csv(path_to_predictDB_predictions, sep="\t", index_col=0)
predictDB_predictions_df.drop("IID", axis=1, inplace=True)
predictDB_predictions_df = predictDB_predictions_df.T

# Load the ground truth counts

count_table_file = parameters["generate_enpact_training_data"]["input_files"]["normalized_expression_data"]
count_df = pd.read_csv(count_table_file, sep="\t", index_col=0)


# Filter for samples, features with missing values

filtered_enpact_df = mean_preds_df.dropna(how="any", axis=0)
filtered_enpact_df = filtered_enpact_df.dropna(how="any", axis=1)
filtered_enpact_df.to_csv(os.path.join(pp_analysis_dir, "filtered_enpact_predictions.txt"), sep="\t")

filtered_predictDB_df = predictDB_predictions_df.dropna(how="any", axis=0)
filtered_predictDB_df.to_csv(os.path.join(pp_analysis_dir, "filtered_predictDB_predictions.txt"), sep="\t")

filtered_ground_truth_df = count_df.dropna(how="any", axis=0)
filtered_ground_truth_df.to_csv(os.path.join(pp_analysis_dir, "filtered_ground_truth.txt"), sep="\t")

# Compute common feature, sample overlap between two predictions and ground truth

non_pd_common_features = list((set(filtered_enpact_df.index) - set(filtered_predictDB_df.index)).intersection(set(filtered_ground_truth_df.index)))

common_features = list(set(filtered_enpact_df.index).intersection(set(filtered_predictDB_df.index).intersection(set(filtered_ground_truth_df.index))))

common_samples = list(set(filtered_enpact_df.columns).intersection(set(filtered_predictDB_df.columns).intersection(set(filtered_ground_truth_df.columns))))

print("Num common features", len(common_features))

print("Num non-predictDB common features", len(non_pd_common_features))

print("Num common samples", len(common_samples))

# Load CV correlations for PredictDB

predictDB_cv_correlations = pd.read_csv(os.path.join(intermediates_dir_predictDB_standard, "predictDB", "database","Model_summary.txt"), sep="\t", index_col=0)

# Subset the dataframes to common features, samples

filtered_enpact_df_matched = filtered_enpact_df.loc[common_features, common_samples]
filtered_enpact_df_nonpd_matched = filtered_enpact_df.loc[non_pd_common_features, common_samples]

filtered_predictDB_df_matched = filtered_predictDB_df.loc[common_features, common_samples]
filtered_cv_correlations = predictDB_cv_correlations.loc[common_features]
filtered_cv_correlations.to_csv(os.path.join(pp_analysis_dir, "filtered_predictDB_cv_correlations.txt"), sep="\t")

filtered_ground_truth_df_matched = filtered_ground_truth_df.loc[common_features, common_samples]
filtered_ground_truth_df_nonpd_matched = filtered_ground_truth_df.loc[non_pd_common_features, common_samples]

# Compute correlation between EnPACT, PredictDB predictions and ground truth

correlation_sets = {
    "EnPACT": {
        "pearson":filtered_enpact_df_matched.corrwith(filtered_ground_truth_df_matched, axis=1),
        "spearman": filtered_enpact_df_matched.corrwith(filtered_ground_truth_df_matched, axis=1, method="spearman"),
        "emp_null_pearson": plotting_utils.empirical_null_correlation(filtered_enpact_df_matched, filtered_ground_truth_df_matched),
        "emp_null_spearman": plotting_utils.empirical_null_correlation(filtered_enpact_df_matched, filtered_ground_truth_df_matched, corr_type="spearman"),
        "anal_null_pearson": plotting_utils.analytical_null_correlation(filtered_enpact_df_matched)

    },
    "EnPACT_nonpd": {
        "pearson": filtered_enpact_df_nonpd_matched.corrwith(filtered_ground_truth_df_nonpd_matched, axis=1),
        "spearman": filtered_enpact_df_nonpd_matched.corrwith(filtered_ground_truth_df_nonpd_matched, axis=1, method="spearman"),
        "emp_null_pearson": plotting_utils.empirical_null_correlation(filtered_enpact_df_nonpd_matched, filtered_ground_truth_df_nonpd_matched),
        "emp_null_spearman": plotting_utils.empirical_null_correlation(filtered_enpact_df_nonpd_matched, filtered_ground_truth_df_nonpd_matched, corr_type="spearman"),
        "anal_null_pearson": plotting_utils.analytical_null_correlation(filtered_enpact_df_nonpd_matched)
    },
    "PredictDB": {
        "pearson": filtered_predictDB_df_matched.corrwith(filtered_ground_truth_df_matched, axis=1),
        "spearman": filtered_predictDB_df_matched.corrwith(filtered_ground_truth_df_matched, axis=1, method="spearman")
    },
    "PredictDB_cv": {
        "pearson": list(filtered_cv_correlations["rho_avg"]),
    },
    "PredictDB_cv_all": {
        "pearson": list(predictDB_cv_correlations["rho_avg"]),
    }

}

correlation_set_pkl_path = os.path.join(pp_analysis_dir, "correlation_sets.pkl")
with open(correlation_set_pkl_path, "wb") as f:
    pkl.dump(correlation_sets, f)

path_to_color_temp = parameters["general_parameters"]["color_palette"]
with open(path_to_color_temp, "r") as f:
    json_color = json.load(f)

correlation_set_colors = {
    "EnPACT": json_color["personalized_training_evaluation"]["EnPACT"],
    "EnPACT_nonpd": json_color["personalized_training_evaluation"]["EnPACT"],
    "PredictDB": json_color["personalized_training_evaluation"]["PredictDB"],
    "PredictDB_cv": json_color["personalized_training_evaluation"]["PredictDB"],
    "PredictDB_cv_all": json_color["personalized_training_evaluation"]["PredictDB"],
    "Null": json_color["personalized_training_evaluation"]["Null"]
}

#############
# START PLOTTING
#############

# Plot cross-individual ground-truth correlations against each other per gene between EnPACT and PredictDB

plot_sets = [
    ["EnPACT", "PredictDB"],
    ["EnPACT", "PredictDB_cv"]
]

for corr_type in ["pearson", "spearman"]:
    for plot_set in plot_sets:

        if corr_type not in correlation_sets[plot_set[0]]:
            continue
        if corr_type not in correlation_sets[plot_set[1]]:
            continue

        corrs_for_plot = pd.DataFrame({plot_set[0]: correlation_sets[plot_set[0]][corr_type], 
                                       plot_set[1]: correlation_sets[plot_set[1]][corr_type]})
        
        sns.scatterplot(data=corrs_for_plot, x=plot_set[0], y=plot_set[1])
        plt.title(f"{plot_set[0]} vs {plot_set[1]} {corr_type} correlation plot")
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.savefig(os.path.join(pp_analysis_plot_dir, f"{plot_set[0]}_{plot_set[1]}_{corr_type}_correlation_plot.png"))
        plt.clf()

        sns.kdeplot(data=corrs_for_plot, x=plot_set[0], y=plot_set[1], fill=True)
        plt.title(f"{plot_set[0]} vs {plot_set[1]} {corr_type} correlation plot")
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.savefig(os.path.join(pp_analysis_plot_dir, f"{plot_set[0]}_{plot_set[1]}_{corr_type}_correlation_plot_kde.png"))
        plt.clf()

        corrs_for_plot = pd.DataFrame({plot_set[0]: np.abs(correlation_sets[plot_set[0]][corr_type]), 
                                       plot_set[1]: correlation_sets[plot_set[1]][corr_type]})
        
        sns.scatterplot(data=corrs_for_plot, x=plot_set[0], y=plot_set[1])
        plt.title(f"{plot_set[0]} vs {plot_set[1]} {corr_type} correlation plot")
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.savefig(os.path.join(pp_analysis_plot_dir, f"{plot_set[0]}_{plot_set[1]}_{corr_type}_abscorrelation_plot.png"))
        plt.clf()

        sns.kdeplot(data=corrs_for_plot, x=plot_set[0], y=plot_set[1], fill=True)
        plt.title(f"{plot_set[0]} vs {plot_set[1]} {corr_type} correlation plot")
        plt.axhline(y=0, color='k')
        plt.axvline(x=0, color='k')
        plt.savefig(os.path.join(pp_analysis_plot_dir, f"{plot_set[0]}_{plot_set[1]}_{corr_type}_abscorrelation_plot_kde.png"))
        plt.clf()


# Plot histograms of the correlation distributions

plot_sets = [
    ["PredictDB_cv_all"],
    ["EnPACT_nonpd"],
    ["EnPACT", "PredictDB_cv"]
]

for plot_set in plot_sets:
    for_hist_plot_dict = {}
    for val in plot_set:
        for_hist_plot_dict[val] = correlation_sets[val]["pearson"]

    for_hist_plot_df = pd.DataFrame(for_hist_plot_dict)
    for_hist_plot_df_flat = for_hist_plot_df.melt()
    for_hist_plot_df_flat["colors"] = [correlation_set_colors[x] for x in for_hist_plot_df_flat["variable"]]
    cur_title = " ".join(plot_set)+" pearsonR histogram"
    cur_savefile = os.path.join(pp_analysis_plot_dir, "_".join(plot_set)+"_pearsonR_hist.png")

    plotting_utils.hist_pearson_R_plot(for_hist_plot_df, 
                                       len(common_samples),
                                        color_mapping=correlation_set_colors, 
                                       title=cur_title, 
                                       path_to_save=cur_savefile)     

# Plot qqplots of the correlations

plot_sets = [
    "EnPACT",
    "EnPACT_nonpd",
    "PredictDB_cv",
    "PredictDB_cv_all"
]


for plot_set in plot_sets:

    cur_title = f"QQPlot of {plot_set}"
    cur_save_file_root = os.path.join(pp_analysis_plot_dir, f"{plot_set}_qqplot")
    plotting_utils.qqR2_python(correlation_sets[plot_set]["pearson"], 
    len(common_samples), 
    title=cur_title,
    path_to_save=cur_save_file_root,
    color=correlation_set_colors[plot_set])


# Plot a few features corrplots

num_features_to_plot = 5

plotted_features = 0
for feature in common_features:

    if plotted_features == num_features_to_plot:
        break

    feature_enpact = list(filtered_enpact_df_matched.loc[feature])
    feature_predictDB = list(filtered_predictDB_df_matched.loc[feature])
    feature_ground_truth = list(filtered_ground_truth_df_matched.loc[feature])

    feature_df = pd.DataFrame({"EnPACT": feature_enpact, "PredictDB": feature_predictDB, "Ground_truth": feature_ground_truth})

    sns.scatterplot(data=feature_df, x="EnPACT", y="Ground_truth", color=correlation_set_colors["EnPACT"])
    plt.title(f"{feature} EnPACT correlation plot")
    plt.savefig(os.path.join(pp_analysis_plot_dir, f"{feature}_EnPACT_correlation_plot.png"))
    plt.clf()

    sns.scatterplot(data=feature_df, x="PredictDB", y="Ground_truth", color=correlation_set_colors["PredictDB"])
    plt.title(f"{feature} PredictDB correlation plot")
    plt.savefig(os.path.join(pp_analysis_plot_dir, f"{feature}_PredictDB_correlation_plot.png"))
    plt.clf()

    plotted_features += 1

