library(tidyverse)
library(httpgd)
library(glmnet)
library(glue)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(stringr)
library(jsonlite)


#################################################################
# 0.) Define input parameters
#################################################################


suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--individual", help='Directory with training data'),
    make_option("--intermediates_dir", help='Directory with intermediate files'),
    make_option("--trained_model", help='.rds file of trained elnet')
)
opt <- parse_args(OptionParser(option_list=option_list))


#################################################################
# 1.) Predict EnPACT prediction for given individual
#################################################################

print(glue("Predicting EnPACT predictions for individual {opt$individual}"))

ind = opt$individual

intermediates_dir <- opt$intermediates_dir

glmnet_model_path <- opt$trained_model
glmnet_model <- readRDS(glmnet_model_path)

linear_dir <- opt$output_dir

X_path <- glue("{intermediates_dir}/personalized_epigenome_predictions_mean_{ind}.txt")
print(ind)
if (file.exists(X_path)) {

    X <- as.matrix(data.table::fread(X_path))

    glmnet_predictions <- as.data.frame(predict(glmnet_model, X, s = "lambda.min", type="response"))

    write_delim(glmnet_predictions, glue("{intermediates_dir}/enpact_personalized_predictions_{ind}.txt"), delim = "\t")

} else {
    print(glue({"{ind} missing"}))
}
