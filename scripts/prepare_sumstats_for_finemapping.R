# Credit to Temi Adeluwa and Sofia Salazar for building the following code for extracting fine-mapped loci

library(optparse)
library(tidyverse)
library(qqman)

#################################################################
# 0.) Define input parameters
#################################################################

option_list <- list(
    make_option("--gwas_sum_stats", help='Path to processed gwas sum stats'),
    make_option("--gwas_sum_stats_name", help='Name of the gwas sum stats'),
    make_option("--output_dir", help='Path to output directory')

)
opt <- parse_args(OptionParser(option_list=option_list))

#################################################################
# 1.) Explore and reformat summary statistics for SusieR
#################################################################

allSumstats <- data.table::fread(opt$gwas_sum_stats)

# remove nonSNPs
cleanSumstats <- allSumstats[nchar(allSumstats$effect_allele) == 1 & nchar(allSumstats$non_effect_allele) == 1, ]

dim(allSumstats)
dim(cleanSumstats)

rm(allSumstats)
# add gwas_significant column

cleanSumstats$gwas_significant <- ifelse(cleanSumstats$pvalue < 5e-08, "YES","NO")

table(cleanSumstats$gwas_significant)

cleanSumstats$zscore <- cleanSumstats$beta / cleanSumstats$standard_error

cleanSumstats$chr_numeric <- as.numeric(str_replace_all(cleanSumstats$chromosome,"chr",""))
table(cleanSumstats$chr_numeric)



colnames(cleanSumstats) <- c("chrom", "SNP", "ref", "alt", "beta", "P", "BP", "frequency", "se", "sample_size", "gwas_significant", "zscore", "CHR")


cleanSumstatsPlotPath <- file.path(opt$output_dir, paste0(opt$gwas_sum_stats_name, ".cleaned.susie.png"))
png(cleanSumstatsPlotPath)
manhattan(na.omit(cleanSumstats))
dev.off()

colnames(cleanSumstats) <- c("chrom", "rsid", "ref", "alt", "beta", "pval", "pos", "frequency", "se", "sample_size", "gwas_significant", "zscore", "chr_numeric")


cleanSumstatsOutputPath <- file.path(opt$output_dir, paste0(opt$gwas_sum_stats_name, ".cleaned.susie.txt"))
write.table(cleanSumstats, file = cleanSumstatsOutputPath, sep = "\t", row.names = F,
            col.names = T,  quote = F)


#################################################################
# 2.) Split sumstats by chromosome
#################################################################

sumstats <- data.table::fread(cleanSumstatsOutputPath)

sumstats$chrom <- gsub("chr", "", sumstats$chrom)
sumstats$chrom <- as.numeric(as.character(sumstats$chrom))
nas <- sumstats[is.na(sumstats$chrom), ]
dim(nas[nas$gwas_significant == "YES"]) # 0 13

filtered_sumstats <- sumstats[!is.na(sumstats$chrom), ]
unique_chroms <- unique(filtered_sumstats$chrom)
for (c in unique_chroms){
  print(c)
  chrom_to_compare <- as.numeric(c)
  chr_filt <- filtered_sumstats[filtered_sumstats$chrom == chrom_to_compare,]
  print(dim(chr_filt))
  write.table(chr_filt, file = paste0(opt$output_dir, "/chr", c, "_snps.txt"), sep = '\t', row.names = F, quote = F)
}

#################################################################
# 3.) Run SusieR
#################################################################
