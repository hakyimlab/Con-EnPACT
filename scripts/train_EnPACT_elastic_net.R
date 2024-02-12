# Author: Temi
# Date: Thursday July 27 2023
# Description: script to train elastic net TFPred models
# Usage: Rscript train_enet.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--train_data_file_X", help='data to train with enet, inputs'),
    make_option("--train_data_file_y", help='data to train with enet, outputs'),
    make_option("--rds_file", help='.rds file to be created as the model'),
    make_option("--nfolds", type="integer", default=5L, help='How many cv folds?'))

opt <- parse_args(OptionParser(option_list=option_list))


library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)
library(tidyverse)

seed <- 2023
if(file.exists(opt$train_data_file_X)){
    print(glue('INFO - Reading train data X...'))
    X_train <- as.matrix(read.table(opt$train_data_file_X, row.names=1, sep='\t', header=F))
    X_train <- X_train[order(row.names(X_train)),]
} else {
    stop(glue('ERROR - Training data cannot be found.'))
}
if(file.exists(opt$train_data_file_y)){
    print(glue('INFO - Reading train data y...'))
    y_train <- as.matrix(read.table(opt$train_data_file_y, row.names=1, sep='\t', header=F))
    y_train <- y_train[order(row.names(y_train)),]
} else {
    stop(glue('ERROR - Training data cannot be found.'))
}

print(glue('INFO - Training data X has {dim(X_train)[1]} rows and {dim(X_train)[2]} columns'))
print(glue('INFO - Training data y has {dim(y_train)[1]} rows and {dim(y_train)[2]} columns'))


cl <- 12 #parallel::makeCluster(5)
print(glue('INFO - Found {parallel::detectCores()} cores but using {cl}'))

set.seed(seed)

doParallel::registerDoParallel(cl)
print(glue('INFO - training enet model'))

# train_methods <- c('linear', 'logistic')
train_methods <- c('linear')

parallel::mclapply(train_methods, function(each_method){

    cl <- 6 #parallel::makeCluster(5)
    doParallel::registerDoParallel(cl)

    print(glue('INFO - Starting to build {each_method} enet model'))

    if(each_method == 'linear'){

        cv_model <- tryCatch({
            glmnet::cv.glmnet(x=X_train, y=y_train, family = "gaussian", type.measure = "mse", alpha = 0.5, keep=T, parallel=T, nfolds=opt$nfolds)
        }, error = function(e){
            print(glue('ERROR - {e}'))
            return(NULL)
        })
        save_name <- paste0(opt$rds_file, '.linear.rds', sep='') #gsub('.rds', '.linear.rds', opt$rds_file)
    } else if (each_method == 'logistic'){

        cv_model <- tryCatch({
            glmnet::cv.glmnet(x=X_train, y=y_train, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=opt$nfolds, trace.it=F)
        }, error = function(e){
            print(glue('ERROR - {e}'))
            return(NULL)
        })
        save_name <-  paste0(opt$rds_file, '.logistic.rds', sep='') #gsub('.rds', '.logistic.rds', opt$rds_file)
    }
    print(cv_model)
    print(glue('INFO - Saving `{save_name}`'))
    #rds_file <- glue('{model_file_basename}.{each_method}.rds')
    saveRDS(cv_model, file=save_name)
    doParallel::stopImplicitCluster()

}, mc.cores=2)

print(glue('INFO - Finished with model training and saving'))
doParallel::stopImplicitCluster()