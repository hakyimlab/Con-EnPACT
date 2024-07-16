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
    make_option("--input_file", help='Path to input file'),
    make_option("--output_file", help='Path of desired output file'),
    make_option("--trained_model", help='.rds file of trained elnet')
)
opt <- parse_args(OptionParser(option_list=option_list))


#################################################################
# 1.) Predict EnPACT prediction for given individual
#################################################################

print(glue("Predicting EnPACT predictions for {opt$input_file}"))

glmnet_model_path <- opt$trained_model
glmnet_model <- readRDS(glmnet_model_path)

if (file.exists(opt$input_file)) {

    X <- as.matrix(data.table::fread(opt$input_file))
    print(dim(X))
    
    X_without_rownames <- X[, -1]
    X_rownames <- X[, 1]

    glmnet_predictions <- as.data.frame(predict(glmnet_model, X_without_rownames, s = "lambda.min", type="response"))
    write_delim(glmnet_predictions, opt$output_file, delim = "\t", col_names = FALSE)
    glmnet_predictions <- cbind(X_rownames, glmnet_predictions)
    write_delim(as.data.frame(X_rownames), paste0(opt$output_file, ".rownames"), delim = "\t", col_names = FALSE)


} else {
    print(glue({"{opt$input_file} missing"}))
}
