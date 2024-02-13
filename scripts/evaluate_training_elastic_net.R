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
    make_option("--train_data_dir", help='Directory with training data'),
    make_option("--rds_file", help='.rds file to be created as the model'),
    make_option("--gene_annotations", help='Annotations file for genes'),
    make_option("--context", help='Context of enpact model'),
    make_option("--epigenome_reference_track", help='Index of reference epigenome track closes match to condition of interest'),
    make_option("--output_dir", help='Output directory to save results'),
    make_option("--color_palette", help='Color palette json file')
)
opt <- parse_args(OptionParser(option_list=option_list))

color_palette <- read_json(opt$color_palette)

#################################################################
# 1.) Define evaluation functions
#################################################################


collect_all_records <- function(model_path, data_subset, 
                                gene_annotations, X_path, y_path, 
                                epi_ref_track) {
  print(glue("********* Collecting data records for data subset: {data_subset}"))

  X <- as.matrix(read.table(X_path, row.names = 1))
  X <- X[order(row.names(X)),]

  print(y_path)
  y <- as.data.frame(read.table(y_path, row.names = 1))
  y$temp <- 1
  y <- y[rownames(X),]
  y$temp <- NULL

  loaded_model <- readRDS(model_path)
  model_predictions <- as.numeric(predict(loaded_model, X, s = "lambda.min", type="response"))

  cur_annotations <- gene_annotations[rownames(y),]

  return_df <- data.frame(data_subset= data_subset, 
    model_predictions=model_predictions, 
    ground_truth = as.numeric(y[,1]),
    epi_ref_track = as.numeric(X[,epi_ref_track]), 
    ensembl_id = cur_annotations[["ensembl_gene_id"]], 
    external_gene_name = cur_annotations[["external_gene_name"]], 
    region = cur_annotations[["target_regions"]], 
    chroms = cur_annotations[["chromosome_name"]])

  return(return_df) 
}


#################################################################
# 2.) Load gene annotations
#################################################################

print("Loading gene annotations")

gene_annotation_file <- opt$gene_annotations
gene_annotation_table <- read.delim(gene_annotation_file, header = TRUE, sep=",")
rownames(gene_annotation_table) <- gene_annotation_table$ensembl_gene_id


#################################################################
# 3.) Run predictions on all data subsets and collect outputs
#################################################################

print("Collecting data records for all data subsets")

all_records_df_store <- list()
c=1

for (data_subset in c("train","test","valid")){
  X_path <- glue(opt$train_data_dir, "/epigenome_{data_subset}.txt")
  y_path <- glue(opt$train_data_dir, "/{data_subset}_{opt$context}_mean_expression.tsv")
  all_records_df_store[[c]] <- collect_all_records(opt$rds_file,
                                                    data_subset, gene_annotation_table, 
                                                    X_path, y_path, as.numeric(opt$epigenome_reference_track))
  c=c+1
}

full_df <- do.call(rbind, all_records_df_store)

full_df <- full_df[!rowSums(is.na(full_df)) > 0,]
full_df$data_subset <- factor(full_df$data_subset, levels = c("train", "valid", "test"))
full_df$chroms <- factor(full_df$chroms, levels = c("all", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"))


full_df$log2_epi_ref_track <- log2(full_df$epi_ref_track)

write.table(full_df, file = glue("{opt$output_dir}/full_predictions_df.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

#################################################################
# 4.) Plot gene counts in each data subset  
#################################################################

print("Plotting gene counts in each data subset")

out_folder <- opt$output_dir

gene_counts_df <- as.data.frame(as.matrix(table(full_df$data_subset)))
gene_counts_df["all","V1"] <- nrow(full_df)
gene_counts_df$data_subset <- rownames(gene_counts_df)
colnames(gene_counts_df) <- c("gene_count", "data_subset")
gene_counts_df$data_subset <- factor(gene_counts_df$data_subset, levels = c("train", "valid", "test", "all"))

png(glue("{out_folder}/counts_data_subsets.png"), width=800, height=400)
ggplot(gene_counts_df) + 
geom_bar(aes(x=data_subset, y=gene_count, fill=data_subset), stat="identity", position="dodge") + 
xlab("Data subset or chromosome") + 
ylab("Gene count") + 
scale_fill_manual(values = unlist(color_palette[["data_subset_colors"]])) + 
theme_minimal(base_size=15)
dev.off()

gene_counts_chrom_df <- distinct(full_df[,c("data_subset", "chroms")])
gene_counts_chrom_df$gene_count <- as.data.frame(as.matrix(table(full_df$chroms)))[gene_counts_chrom_df$chroms,]

png(glue("{out_folder}/counts_all_chroms_both.png"), width=800, height=400)
ggplot(gene_counts_chrom_df) + 
  geom_bar(aes(x=chroms, y=gene_count, fill=data_subset), stat="identity", position="dodge") + 
  xlab("Data subset or chromosome") + 
  ylab("Gene count") + 
  facet_wrap(~data_subset, ncol=1) + 
  scale_fill_manual(values = unlist(color_palette[["data_subset_colors"]]))+ 
  theme_minimal(base_size=15)
dev.off()


#################################################################
# 4.) Plot predictive performance across genes  
#################################################################

print("Plotting predictive performance across genes")

data_subsets <- c("train", "valid", "test")
cor_names <- c("EnPACT_pred_vs_ground_truth", "log2_epi_ref_track_vs_ground_truth")
cor_vecs <- c("model_predictions", "log2_epi_ref_track")
cor_against <- c("ground_truth", "ground_truth")
cor_types <- c("spearman","pearson")

corrs_df_for_subsets <- matrix(nrow = length(data_subsets)*length(cor_vecs)*length(cor_types), ncol = 4)

i = 1
for (data_subset in data_subsets) {
  cur_subset_df <- full_df[full_df$data_subset == data_subset,]
  cur_subset_df <- cur_subset_df[!rowSums(is.na(cur_subset_df)) > 0,]
  for (cor_type in cor_types) {
    for (c in 1:length(cor_names)) {
      cur_cor <- cor(cur_subset_df[,cor_vecs[c]], cur_subset_df[,cor_against[c]], method = cor_type)
      corrs_df_for_subsets[i,] <- c(data_subset, cor_names[c], cor_type, cur_cor)
      i = i + 1
    }
  }
}

corrs_df_for_subsets <- as.data.frame(corrs_df_for_subsets)
colnames(corrs_df_for_subsets) <- c("data_subset", "cor_id", "cor_type", "corrs")
corrs_df_for_subsets$cor_id <- factor(corrs_df_for_subsets$cor_id, levels = cor_names)
corrs_df_for_subsets$corrs <- as.numeric(as.character(corrs_df_for_subsets$corrs))
corrs_df_for_subsets$data_subset <- factor(corrs_df_for_subsets$data_subset, levels = c("train", "valid", "test"))

png(glue("{out_folder}/all_corrs_subsets.png"), width=600, height=500)
ggplot(corrs_df_for_subsets) + 
  geom_bar(aes(x=data_subset, y=corrs, fill=cor_id), stat="identity", position="dodge") + 
  xlab("Data subset") + ylab("Correlation") + 
  facet_wrap(~cor_type, ncol=1) + 
  scale_fill_manual(values = unlist(color_palette[["enpact_training_evaluation"]]))+ 
  theme_minimal(base_size=15)
dev.off()
