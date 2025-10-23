## This script generates figure 2 and associated supplementary files - plots for different substitution types, distribution of editing sites across different genic regions, number of exonic substitution types, Alu editing, correlation plots, plots of unique edits, overrepresentation plots for GO-BP,Reactome, and AEI plot.

# Figure 2
#libraries
library(dplyr)
library(stringr)
library(tidyverse)

#============  copy files =======================================================

# copy the .vcf files from the VCf folder to a new folder for analysis
source_dir <- "G:\\Manuscript_2_rerunning_new_SG\\Original_Data\\vcf_files\\vcf_files\\final"
destination_dir <- "G:\\Manuscript_2_rerunning_new_SG\\Total_substitutions"
files <- list.files(source_dir, pattern = "SRR.*filtered_snps_finalAllNoSnpsediting\\.vcf", full.names = TRUE)
file.copy(files, destination_dir, overwrite = TRUE)

#=========== read vcf by skipping lines starting with ## and save as csv =======

# Directory containing VCF files
vcf_dir <- "G:\\Manuscript_2_rerunning_new_SG\\Total_substitutions"  

# List of all  VCF files in the above folder
vcf_files <- list.files(vcf_dir, pattern = "\\.vcf$", full.names = TRUE)

# This will Loop through each VCF file to read files, remove ## lines, extract header from the first remaining line, and read rest of the data 
# outpu saved as csv

for (vcf_file in vcf_files) {
  vcf <- readLines(vcf_file)
  vcf <- vcf[!grepl("^##", vcf)]
  header <- strsplit(vcf[1], "\t")[[1]]
  data <- read.table(text = paste(vcf[-1], collapse = "\n"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(data) <- header
  output_csv <- sub("\\.vcf$", ".csv", vcf_file)
  write.csv(data, output_csv, row.names = FALSE)
  print(paste("Processed:", vcf_file, "Saved:", output_csv))
}

#============  copy file =======================================================

# copy the .csv files from the VCf folder to a new folder for analysis
source_dir <- "G:\\Manuscript_2_rerunning_new_SG\\Total_substitutions"
destination_dir <- "G:\\Manuscript_2_rerunning_new_SG\\Total_substitutions\\csv\\"
files <- list.files(source_dir, pattern = "SRR.*filtered_snps_finalAllNoSnpsediting\\.csv", full.names = TRUE)
file.copy(files, destination_dir, overwrite = TRUE)

#=========== Applying filters ============================================================
# Working directory
setwd("G:\\Manuscript_2_rerunning_new_SG\\Total_substitutions\\csv\\")

# This function will run on each file "*filtered_snps_finalAll.csv" and apply the filters
process_file <- function(file_path) {
  
  # Adding column names
  col_names <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "20", "k", "LINK", "Substitutions")
  VCF_files <- read.csv(file_path, header = FALSE, sep = ",", col.names = col_names)
  
  # Add a LINK column (to use later)
  VCF_files <- VCF_files %>% mutate(LINK = "-to-")
  
  # Merge the three columns "ref" link" "alt"  to make different substitutio types
  VCF_files <- VCF_files %>% mutate(Substitutions = paste(REF,LINK, ALT, sep = ""))
  
  # Filter1: keep rows with "." in "ID" column, this is to remove the rows with SNP ids
  filtered_data <- VCF_files %>%
    filter(ID == ".")
  
  # Filter2: remove rows that have 2 values in column "ALT" - ie these are polymorphic sites
  filtered_data <- filtered_data %>%
    filter(!str_detect(ALT, ","))
  
  # Filter3: keep rows where either "REF = A" and "ALT = G" or "REF = T" and "ALT = C" - to keep only ADAR edited sites
  #filtered_data <- filtered_data %>%
  #filter((REF == "A" & ALT == "G") | (REF == "T" & ALT == "C"))
  
  # Count total number each substitutions
  ref_tally_AC <- sum(filtered_data$Substitutions == "A-to-C")
  ref_tally_AG <- sum(filtered_data$Substitutions == "A-to-G")
  ref_tally_AT <- sum(filtered_data$Substitutions == "A-to-T")
  ref_tally_CA <- sum(filtered_data$Substitutions == "C-to-A")
  ref_tally_CG <- sum(filtered_data$Substitutions == "C-to-G")
  ref_tally_CT <- sum(filtered_data$Substitutions == "C-to-T")
  ref_tally_GA <- sum(filtered_data$Substitutions == "G-to-A")
  ref_tally_GC <- sum(filtered_data$Substitutions == "G-to-C")
  ref_tally_GT <- sum(filtered_data$Substitutions == "G-to-T")
  ref_tally_TA <- sum(filtered_data$Substitutions == "T-to-A")
  ref_tally_TC <- sum(filtered_data$Substitutions == "T-to-C")
  ref_tally_TG <- sum(filtered_data$Substitutions == "T-to-G")
  #ref_tally_T <- sum(filtered_data$REF == "T") # will include all the T to C
  
  # Extract the file name without extension to identify files in the final dataframe and save
  file_name <- tools::file_path_sans_ext(basename(file_path)) ## tools package
  
  write.csv(VCF_files, file_path, row.names = FALSE)
  
  return(list(ref_tally_AC = ref_tally_AC, ref_tally_AG = ref_tally_AG, ref_tally_AT = ref_tally_AT,
              ref_tally_CA = ref_tally_CA, ref_tally_CG = ref_tally_CG, ref_tally_CT = ref_tally_CT,
              ref_tally_GA = ref_tally_GA, ref_tally_GC = ref_tally_GC, ref_tally_GT = ref_tally_GT, 
              ref_tally_TA = ref_tally_TA, ref_tally_TC = ref_tally_TC, ref_tally_TG = ref_tally_TG,file_name = file_name))
}

# Directory containing the files
directory_path <- "G:\\Manuscript_2_rerunning_new_SG\\Total_substitutions\\csv\\"

# Empty list for results 
results_list <- list()

# Loop through the files
files <- list.files(pattern = "*.csv", full.names = TRUE, path = directory_path)
for (file_path in files) {
  result <- process_file(file_path)
  results_list[[length(results_list) + 1]] <- result # to ensure that the result is added to the end of the list
}

# Extracting vectors of A, T, and file names 
ref_tally_AC_vector <- sapply(results_list, function(result) result$ref_tally_AC)
ref_tally_AG_vector <- sapply(results_list, function(result) result$ref_tally_AG)
ref_tally_AT_vector <- sapply(results_list, function(result) result$ref_tally_AT)
ref_tally_CA_vector <- sapply(results_list, function(result) result$ref_tally_CA)
ref_tally_CG_vector <- sapply(results_list, function(result) result$ref_tally_CG)
ref_tally_CT_vector <- sapply(results_list, function(result) result$ref_tally_CT)
ref_tally_GA_vector <- sapply(results_list, function(result) result$ref_tally_GA)
ref_tally_GC_vector <- sapply(results_list, function(result) result$ref_tally_GC)
ref_tally_GT_vector <- sapply(results_list, function(result) result$ref_tally_GT)
ref_tally_TA_vector <- sapply(results_list, function(result) result$ref_tally_TA)
ref_tally_TC_vector <- sapply(results_list, function(result) result$ref_tally_TC)
ref_tally_TG_vector <- sapply(results_list, function(result) result$ref_tally_TG)
file_names_vector <- sapply(results_list, function(result) result$file_name)

# Joining all the above vectors of number of diffrent substitutions to a dataframe
Total_substitutions <- data.frame(sample = file_names_vector, A_to_C = ref_tally_AC_vector, A_to_G = ref_tally_AG_vector, 
                                  A_to_T = ref_tally_AT_vector,C_to_A = ref_tally_CA_vector,
                                  C_to_G = ref_tally_CG_vector,C_to_T = ref_tally_CT_vector,
                                  G_to_A = ref_tally_GA_vector,G_to_C = ref_tally_GC_vector,
                                  G_to_T = ref_tally_GT_vector,T_to_A = ref_tally_TA_vector,
                                  T_to_C = ref_tally_TC_vector,T_to_G = ref_tally_TG_vector)
# Add condition and write
Total_substitutions$Condition <- c(rep("Non-critical", 23), rep("Critical", 46))
write.csv(Total_substitutions, "Result_Total_substitutions.csv", row.names = FALSE)

#============ PLOT  (Fig 2A)=============================================================

#libraries
library(ggplot2)
library(reshape2)

# Data
Substitution_data <- read.csv("Result_Total_substitutions.csv", header = TRUE)

# Convert this data to long format
Substitution_longFormat <- melt(Substitution_data, id.vars = c("sample", "Condition"), 
                                measure.vars = c("A_to_G", "T_to_C", "A_to_C","A_to_T", 
                                                 "C_to_A", "C_to_G", "C_to_T", 
                                                 "G_to_A", "G_to_C", "G_to_T", 
                                                 "T_to_A","T_to_G"))

# Plot
plotA <- ggplot(Substitution_longFormat, aes(x = variable, y = value, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = c("tomato", "darkblue")) + # Custom colors
  labs(x = "Type of RNA modification",
       y = "Number of RNA editing sites") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14), 
        legend.title = element_blank(),  # This will remove the title of legend
        legend.position = c(0.9, 0.9),  # This will position legend inside on the top corner
        legend.justification = c("right", "top")) 

ggsave("Total_substitutions_AIDD_VCF_ggsave.png", width = 8, height = 5)

#======================= only ADAR edits=============================================
# Read CSV file
Total_subst <- read.csv("Result_Total_substitutions.csv")  # Adjust separator if needed

# Filter and calculate ADAR total
adar_Total_subst <- Total_subst %>%
  select(sample, A_to_G, T_to_C, Condition) %>%
  mutate(ADAR_total = A_to_G + T_to_C)

# Save only ADAR substitution file
write.csv(adar_Total_subst, "Total_ADAR_only_substitutions.csv", row.names = FALSE)


# Calculate mean values per group
means_adar_Total_subst <- adar_Total_subst %>%
  group_by(Condition) %>%
  summarise(mean_value = mean(ADAR_total))

# Plot with custom colors
ggplot(adar_Total_subst, aes(x = Condition, y = ADAR_total, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +
  geom_text(data = means_adar_Total_subst, aes(x = Condition, y = mean_value, label = round(mean_value, 1)),
            vjust = -0.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Critical" = "tomato", "Non-critical" = "darkblue")) +
  labs(y = "Total number of ADAR substitutions",
       x = "Condition") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14), 
        legend.title = element_blank(),  
        legend.position = ("none"))  
        
ggsave("Total_substitutions_AIDD_VCF_ADAR_only.png", width = 8, height = 5)

#-------------------------------------------------------------------------------
# correlation between ADAR expression and editing 
#-------------------------------------------------------------------------------

library(tidyverse)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Define folder path where all TPM files are stored
ADAR_TPM_counts_dir <- "G:\\Manuscript_2_rerunning_new_SG\\Total_substitutions\\csv\\ADAR_counts" 

# List all files
ADAR_expression_files <- list.files(ADAR_TPM_counts_dir, full.names = TRUE)

# Initialize an empty list to store TPMs
TPM_list <- list()

# Function to extract TPMs for ADAR genes from one file
extract_ADAR_TPMs <- function(file_path) {
  df <- read.csv(file_path, sep = "\t")  
  
  # Extract sample name 
  sample_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Filter for ADAR genes 
  df_filtered <- df %>%
    filter(`Gene.Name` %in% c("ADAR", "ADARB1", "ADARB2")) %>%
    select(`Gene.Name`, TPM)
  
  # Reshape to wide format
  df_wide <- df_filtered %>%
    pivot_wider(names_from = `Gene.Name`, values_from = TPM) %>%
    mutate(sample = sample_name) %>%
    select(sample, ADAR, ADARB1, ADARB2)  # Rename later if needed
}

# Apply function to all files
tpm_list <- lapply(ADAR_expression_files, extract_ADAR_TPMs)

# Combine all into one dataframe
adar_expression_df <- bind_rows(tpm_list)

# Rename for clarity
adar_expression_df <- adar_expression_df %>%
  rename(ADAR1_TPM = ADAR, ADAR2_TPM = ADARB1, ADAR3_TPM = ADARB2)
write.csv(adar_expression_df, "All_adar_expression.csv", row.names = FALSE)

# load data --------------------------------------------------------------------

# combine ADAR exxpression and editing files
## Note: manually remove the extra from the sample names in both the substitution and expression file by keeping only SRR numbers 
# Load both the dataframes (adar expression and total number of ADAR edits)

adar_Total_subst <- read.csv("Total_ADAR_only_substitutions.csv")
adar_expression_df <- read.csv("All_adar_expression.csv")

# Merge expression and editing data by sample names (SRR)
adar_combined <- left_join(adar_Total_subst, adar_expression_df, by = "sample")
# 
# # Log2 transform TPM values
# adar_combined <- adar_combined %>%
#   mutate(log_ADAR1_TPM = log2(ADAR1_TPM + 1),
#          log_ADAR2_TPM = log2(ADAR2_TPM + 1),
#          log_ADAR3_TPM = log2(ADAR3_TPM + 1))


# Spearman correlation by condition to check if there is any association
# adar_combined %>%
#   group_by(Condition) %>%
#   summarise(
#     cor_ADAR1 = cor(ADAR_total, ADAR1_TPM, method = "spearman"),
#     cor_ADAR2 = cor(ADAR_total, ADAR2_TPM, method = "spearman"),
#     cor_ADAR3 = cor(ADAR_total, ADAR3_TPM, method = "spearman")
#   )


# Linear regression ------------------------------------------------------------

# ADAR1
fit_critical_ADAR1 <- lm(ADAR_total ~ ADAR1_TPM, data = subset(adar_combined, Condition == "Critical"))
summary(fit_critical_ADAR1)
fit_noncritical_ADAR1 <- lm(ADAR_total ~ ADAR1_TPM, data = subset(adar_combined, Condition == "Non-critical"))
summary(fit_noncritical_ADAR1)

# ADAR2
fit_critical_ADAR2 <- lm(ADAR_total ~ ADAR2_TPM, data = subset(adar_combined, Condition == "Critical"))
summary(fit_critical_ADAR2)
fit_noncritical_ADAR2 <- lm(ADAR_total ~ ADAR2_TPM, data = subset(adar_combined, Condition == "Non-critical"))
summary(fit_noncritical_ADAR2)

# ADAR3
fit_critical_ADAR3 <- lm(ADAR_total ~ ADAR3_TPM, data = subset(adar_combined, Condition == "Critical"))
summary(fit_critical_ADAR3)
fit_noncritical_ADAR3 <- lm(ADAR_total ~ ADAR3_TPM, data = subset(adar_combined, Condition == "Non-critical"))
summary(fit_noncritical_ADAR3)

# plot ADAR1 --------------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(ggpmisc)

ggplot(adar_combined, aes(x = ADAR1_TPM, y = ADAR_total)) +
  geom_point(aes(color = Condition), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, aes(color = Condition)) +
  
    # Linear model R2 and pval
  stat_fit_glance(
    method = "lm",
    method.args = list(formula = y ~ x),
    aes(label = paste0("adj R² = ", signif(..adj.r.squared.., 3), 
                       ", p = ", signif(..p.value.., 3))),
    hjust = -0.1, vjust = 1.4
  ) +
  
  
  facet_wrap(~Condition) +
  labs(x = "ADAR1 expression (TPM)",
       y = "Number of ADAR edits") +
  scale_color_manual(values = c("Critical" = "tomato", "Non-critical" = "darkblue")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none",
        strip.text = element_text(size = 16))

ggsave("ADAR1corelation.png", width = 8, height = 5)

# plot ADAR2 --------------------------------------------------------------------
ggplot(adar_combined, aes(x = ADAR2_TPM, y = ADAR_total)) +
  geom_point(aes(color = Condition), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, aes(color = Condition)) +
  
    
  # Linear model R2 and pVal
  stat_fit_glance(
    method = "lm",
    method.args = list(formula = y ~ x),
    aes(label = paste0("adj R² = ", signif(..adj.r.squared.., 3), 
                       ", p = ", signif(..p.value.., 3))),
     hjust = -0.1, vjust = 1.4
  ) +
  
 
  facet_wrap(~Condition) +
  labs(x = "ADAR2 expression (TPM)",
       y = "Number of ADAR edits") +
  scale_color_manual(values = c("Critical" = "tomato", "Non-critical" = "darkblue")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none",
        strip.text = element_text(size = 16))
ggsave("ADAR2corelation.png", width = 8, height = 5)

# plot ADAR3 --------------------------------------------------------------------
ggplot(adar_combined, aes(x = ADAR3_TPM, y = ADAR_total)) +
  geom_point(aes(color = Condition), size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, aes(color = Condition)) +
  
  # Linear model R2 and pVal 
  stat_fit_glance(
    method = "lm",
    method.args = list(formula = y ~ x),
    aes(label = paste0("adj R² = ", signif(..adj.r.squared.., 3), 
                       ", p = ", signif(..p.value.., 3))),
    hjust = -0.1, vjust = 1.4
  ) +
  
  
  facet_wrap(~Condition) +
  labs(x = "ADAR3 expression (TPM)",
       y = "Number of ADAR edits") +
  scale_color_manual(values = c("Critical" = "tomato", "Non-critical" = "darkblue")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none",
        strip.text = element_text(size = 16))

ggsave("ADAR3corelation.png", width = 8, height = 5)

#==================== Annovar annotation - genic region ========================

# copy the .csv files from the annovar folder to a new folder for analysis
source_dir <- "G:\\Manuscript_2_rerunning_new_SG\\Original_Data\\annovar_gene_annotation"
destination_dir <- "G:\\Manuscript_2_rerunning_new_SG\\Annovar_Gene_annotation"
files <- list.files(source_dir, pattern = "SRR.*filtered_snps_finalAllNoSnpsediting.hg19_multianno.txt", full.names = TRUE)
file.copy(files, destination_dir, overwrite = TRUE)

#=============== annotate ======================================================

library(dplyr)

setwd("G:\\Manuscript_2_rerunning_new_SG\\Annovar_Gene_annotation") 

# Function to filter and sum up each genic region
process_file <- function(file_path) {
  ADARonly_VCF_files <- read.table(file_path, header = TRUE, sep = "\t")
  
  # Filter1: keep rows where either "REF = A" and "ALT = G" or "REF = T" and "ALT = C" - to keep only ADAR edited sites
  filtered_data <- ADARonly_VCF_files %>% filter((Ref == "A" & Alt == "G") | (Ref == "T" & Alt == "C"))
  
  # Count total number of each substitution
  ref_tally_downstream_vector <- sum(filtered_data$Func.refGene == "downstream", na.rm = TRUE)
  ref_tally_EXON_vector <- sum(filtered_data$Func.refGene == "exonic", na.rm = TRUE)
  ref_tally_intergenic <- sum(filtered_data$Func.refGene == "intergenic", na.rm = TRUE)
  ref_tally_intronic <- sum(filtered_data$Func.refGene == "intronic", na.rm = TRUE)
  ref_tally_SPLICE_SITE_vector <- sum(filtered_data$Func.refGene == "splicing", na.rm = TRUE)
  ref_tally_3UTR <- sum(filtered_data$Func.refGene == "UTR3", na.rm = TRUE)
  ref_tally_5UTR <- sum(filtered_data$Func.refGene == "UTR5", na.rm = TRUE)
  ref_tally_upstream_vector <- sum(filtered_data$Func.refGene == "upstream", na.rm = TRUE)
  
  # Extract file name without extension
  file_name <- tools::file_path_sans_ext(basename(file_path)) ## tools package
  
  return(list(
    ref_tally_intergenic = ref_tally_intergenic, 
    ref_tally_intronic = ref_tally_intronic,
    ref_tally_3UTR = ref_tally_3UTR, 
    ref_tally_5UTR = ref_tally_5UTR,
    ref_tally_EXON_vector = ref_tally_EXON_vector,
    ref_tally_downstream_vector = ref_tally_downstream_vector,
    ref_tally_upstream_vector = ref_tally_upstream_vector,
    ref_tally_SPLICE_SITE_vector = ref_tally_SPLICE_SITE_vector
  ))
}

# Directory containing the vcf files saved as csv
directory_path <- "G:\\Manuscript_2_rerunning_new_SG\\Annovar_Gene_annotation"

# Open an empty list
results_list <- list()

# Loop over each file
files <- list.files(pattern = "*.txt", full.names = TRUE, path = directory_path)
for (file_path in files) {
  result <- process_file(file_path)
  results_list[[length(results_list) + 1]] <- result # Ensures results are added to the list
}

# Extract vectors from results list
ref_tally_intergenic_vector <- sapply(results_list, function(result) result$ref_tally_intergenic)
ref_tally_intronic_vector <- sapply(results_list, function(result) result$ref_tally_intronic)
ref_tally_3UTR_vector <- sapply(results_list, function(result) result$ref_tally_3UTR)
ref_tally_5UTR_vector <- sapply(results_list, function(result) result$ref_tally_5UTR)
ref_tally_SPLICE_SITE_vector <- sapply(results_list, function(result) result$ref_tally_SPLICE_SITE_vector)
ref_tally_EXON_vector <- sapply(results_list, function(result) result$ref_tally_EXON_vector)
ref_tally_downstream_vector <- sapply(results_list, function(result) result$ref_tally_downstream_vector)
ref_tally_upstream_vector <- sapply(results_list, function(result) result$ref_tally_upstream_vector)

# Extract file names
file_names_vector <- sapply(files, function(file) tools::file_path_sans_ext(basename(file)))

# Join all the vectors into a dataframe
Total_substitutions <- data.frame(
  sample = file_names_vector,
  intergenic_region = ref_tally_intergenic_vector, 
  intronic_region = ref_tally_intronic_vector, 
  ref_tally_3UTR = ref_tally_3UTR_vector,
  ref_tally_5UTR = ref_tally_5UTR_vector,
  splice_site = ref_tally_SPLICE_SITE_vector,
  exon = ref_tally_EXON_vector,  
  downstream = ref_tally_downstream_vector,
  upstream = ref_tally_upstream_vector
)


# Add condition and write 
Total_substitutions$Condition <- c(rep("Non-critical", 23), rep("Critical", 46))
write.csv(Total_substitutions, "Result_regionwise_Total_substitutions_with_condition.csv", row.names = FALSE)

#======================================== plot =================================
library(ggplot2)
library(reshape2)
library(dplyr)

# Read the file
region_wise_ADAR_edits <- read.csv("G:\\Manuscript_2_rerunning_new_SG\\Annovar_Gene_annotation\\Result_regionwise_Total_substitutions_with_condition.csv")

# Define region columns
region_columns <- c("intronic_region", "intergenic_region", "ref_tally_3UTR", 
                    "downstream", "exon", "upstream", "ref_tally_5UTR", "splice_site")

# Calculate mean number of edits per region for each condition separately
mean_edits_per_region <- region_wise_ADAR_edits %>%
  group_by(Condition) %>%
  summarise(across(all_of(region_columns), mean, na.rm = TRUE))

# Converting to long format for plotting
long_mean_data <- melt(mean_edits_per_region, id.vars = "Condition", 
                       variable.name = "Region", value.name = "Mean_Edits")

# Rename the regions
region_names <- c("intronic_region" = "Intronic", 
                  "intergenic_region" = "Intergenic", 
                  "ref_tally_3UTR" = "3' UTR", 
                  "downstream" = "Downstream", 
                  "exon" = "Exonic", 
                  "upstream" = "Upstream", 
                  "ref_tally_5UTR" = "5' UTR", 
                  "splice_site" = "Splice Site")

# Applying renaming
long_mean_data$Region <- factor(long_mean_data$Region, levels = names(region_names), 
                                labels = region_names)

# Save the updated dataset
write.csv(long_mean_data, "Mean_regionwise_ADAR_edits_per_condition.csv", row.names = FALSE)

# bar plot
ggplot(long_mean_data, aes(x = Region, y = Mean_Edits, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("tomato", "darkblue")) +  
  labs(y = "Mean Number of ADAR Editing Sites", x = "Genomic Region") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14), 
        legend.title = element_blank(), 
        legend.position = c(0.9, 0.9))  

ggsave("mean_regionwise_ADAR_edits_per_conditions_ggsave.png", width = 8, height = 5)

#================= synonymous and non synonymous exonic ========================

# Read the file
exonic_mutations <- read.csv("G:\\Manuscript_2_rerunning_new_SG\\Annovar_Gene_annotation\\Result_regionwise_Total_substitutions_exonic_syn.csv")

# Mutation type columns
mutation_columns <- c("nonsynonymous", "startloss", "stopgain", "stoploss", "synonymous", "unknown")

# Calculate mean number of mutations per type for each condition separately
mean_mutations_per_type <- exonic_mutations %>%
  group_by(Condition) %>%
  summarise(across(all_of(mutation_columns), mean, na.rm = TRUE))

# Convert to long format for plotting
long_mean_data <- melt(mean_mutations_per_type, id.vars = "Condition", 
                       variable.name = "Mutation_Type", value.name = "Mean_Mutations")

# Rename 
mutation_names <- c("nonsynonymous" = "Nonsynonymous",
                    "synonymous" = "Synonymous",
                    "unknown" = "Unknown",
                    "stoploss" = "Stop Loss",
                    "startloss" = "Start Loss", 
                    "stopgain" = "Stop Gain"
)

# Apply renaming
long_mean_data$Mutation_Type <- factor(long_mean_data$Mutation_Type, levels = names(mutation_names), 
                                       labels = mutation_names)
write.csv(long_mean_data, "Mean_exonic_mutations_per_condition.csv", row.names = FALSE)

# Barplot
ggplot(long_mean_data, aes(x = Mutation_Type, y = Mean_Mutations, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("tomato", "darkblue")) +  
  labs(y = "Mean number of ADAR editing sites", x = "Exonic substitution Type") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 16),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14), 
        legend.title = element_blank(), # Remove legend title
        legend.position = c(0.9, 0.9))  # Position legend on top


ggsave("G:\\Manuscript_2_rerunning_new_SG\\Annovar_Gene_annotation\\Exonic_edit_Distribution_ggsave.png", width = 8, height = 5)

#-------------------------------------------------------------------------------

# Finding unique and shared edits from AIDD vcf annoated with ANNOVAR

#-------------------------------------------------------------------------------

# set WD 
setwd("G:\\Manuscript_2_rerunning_new_SG\\Unique_edits_from_AIDD_VCF_files_annotated_with_Annovar\\")

# libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggvenn)
library(ggplot2)
library(cowplot)

# Critical and NC samples (raw AIDD VCFs annoated with ANNOVAR) are in separate folders
condition_Critical_path <- "G:\\Manuscript_2_rerunning_new_SG\\Unique_edits_from_AIDD_VCF_files_annotated_with_Annovar\\Critical"
condition_Non_critical_path <- "G:\\Manuscript_2_rerunning_new_SG\\Unique_edits_from_AIDD_VCF_files_annotated_with_Annovar\\Non_Critical"

# Function filters ADAR edits and add CHROM_POS column to each file 
process_vcf <- function(file) {
  df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  df <- df %>% filter((Ref == "A" & Alt == "G") | (Ref == "T" & Alt == "C")) 
  df <- df %>% mutate(CHROM_POS = paste(Chr, Start, sep = "_"))
  return(df)
}

# Process all files in a specific condition folder
process_condition <- function(folder_path) {
  files <- list.files(folder_path, full.names = TRUE, pattern = "\\.vcf|\\.txt|\\.tsv|\\.csv")
  all_data <- lapply(files, process_vcf) # apply the process_vcf function to all files in the list 
  names(all_data) <- basename(files)  # getting file names
  return(all_data)
}


# Get the list of CHROM_POS per sample
extract_chrom_pos_lists <- function(all_data) {
  lapply(all_data, function(df) unique(df$CHROM_POS))
}

# Filter CHROM_POS that is present in all the samples ie high confidence
get_high_conf_sites <- function(all_data, chrom_pos_lists) {
  total_samples <- length(chrom_pos_lists)
  shared_sites <- Reduce(intersect, chrom_pos_lists) # "Reduce" and "Intersect" will get sites that are present in all samples
  combined_df <- bind_rows(all_data) %>% as_tibble() # bind all files into one frame
  print(colnames(combined_df))
  print(paste("Shared sites count:", length(shared_sites))) #check
  filtered_df <- combined_df %>% filter(.data$CHROM_POS %in% shared_sites) # Filter the big combined df to include only shared sites
  final_df <- filtered_df %>% distinct(CHROM_POS, .keep_all = TRUE) # Remove duplicates
  return(final_df)
}

# Extract annotations for shared and unique sites
get_annotations_for_sites <- function(df, chrom_pos_list) {
  # Filter the dataframe to include only the rows corresponding to the provided CHROM_POS list
  annotations_df <- df %>% filter(CHROM_POS %in% chrom_pos_list)
  return(annotations_df)
}


# Apply "process_condition" function
data_Critical <- process_condition(condition_Critical_path)
data_Non_critical <- process_condition(condition_Non_critical_path)

# To the above apply "extract_chrom_pos_list" function
chrom_pos_Critical <- extract_chrom_pos_lists(data_Critical)
chrom_pos_Non_critical <- extract_chrom_pos_lists(data_Non_critical)

# Filter high confidence editing sites (present in all samples) - apply "get_high_conf_sites" function
high_conf_Critical <- get_high_conf_sites(data_Critical, chrom_pos_Critical)
high_conf_Non_critical <- get_high_conf_sites(data_Non_critical, chrom_pos_Non_critical)

# Compare: shared & unique
shared_sites <- intersect(high_conf_Critical$CHROM_POS, high_conf_Non_critical$CHROM_POS)
unique_Critical <- setdiff(high_conf_Critical$CHROM_POS, high_conf_Non_critical$CHROM_POS)
unique_Non_critical <- setdiff(high_conf_Non_critical$CHROM_POS, high_conf_Critical$CHROM_POS)

# Get annotations for shared and unique sites
shared_df <- bind_rows(
  high_conf_Critical %>% filter(CHROM_POS %in% shared_sites),
  high_conf_Non_critical %>% filter(CHROM_POS %in% shared_sites)
) %>% distinct(CHROM_POS, .keep_all = TRUE)

dim(shared_df)

unique_annotations_Critical <- bind_rows(
  high_conf_Critical %>% filter(CHROM_POS %in% unique_Critical),
  high_conf_Non_critical %>% filter(CHROM_POS %in% unique_Critical)
)%>% distinct(CHROM_POS, .keep_all = TRUE)
dim(unique_annotations_Critical)

unique_annotations_Non_critical <- bind_rows(
  high_conf_Critical %>% filter(CHROM_POS %in% unique_Non_critical),
  high_conf_Non_critical %>% filter(CHROM_POS %in% unique_Non_critical)
)%>% distinct(CHROM_POS, .keep_all = TRUE)
dim(unique_annotations_Non_critical)


# Save sa CSVs
write.csv(shared_df, "shared_edits_crit_non_crit_present_100.csv", row.names = FALSE)
write.csv(unique_annotations_Critical, "unique_sites_annotated_Critical.csv", row.names = FALSE)
write.csv(unique_annotations_Non_critical, "unique_sites_annotated_Non_critical.csv", row.names = FALSE)

# === Venn Diagram - shared/unique sites========================================
venn_data <- list(
  Critical = high_conf_Critical$CHROM_POS,
  Non_critical = high_conf_Non_critical$CHROM_POS
)
ggvenn(venn_data, fill_color = c("tomato", "darkblue"))
ggsave("shared_and_unique_edits_across_all_samples_2.png", width = 8, height = 7, bg = "white")

# === Venn Diagram - shared/unique genes========================================
venn_data <- list(
  Critical = high_conf_Critical$Gene.refGene,
  Non_critical = high_conf_Non_critical$Gene.refGene
)
ggvenn(venn_data, fill_color = c("tomato", "darkblue"))
ggsave("shared_and_unique_genes_across_all_samples_2.png", width = 8, height = 7, bg = "white")

#========= genomic region_critical_unique=======================================

# Group by genomic region and count the number of unique CHROM_POS 
edit_counts <- unique_annotations_Critical %>%
  group_by(Func.refGene) %>%
  summarize(n_edits = n_distinct(CHROM_POS)) %>%
  arrange(desc(n_edits)) 

# Convert genomic region to a factor 
edit_counts$Func.refGene <- factor(edit_counts$Func.refGene, levels = edit_counts$Func.refGene[order(edit_counts$n_edits, decreasing = TRUE)])

# Plot
Unique_critical_genomic_region <- ggplot(edit_counts, aes(x = Func.refGene, y = n_edits)) +
  geom_bar(stat = "identity", fill = "tomato") +
  geom_text(aes(label = n_edits), vjust = -0.3, size = 3) + # Add labels on top of bars
  labs(
    x = "Genomic Region",
    y = "Number of unique edits - critical"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14),
        legend.position = c(0.9, 0.9),  
        legend.justification = c("right", "top"))

ggsave("Unique_critical_genomic_region.png", width = 8, height = 7)


#---- exonic edits with synonymous and non synonymous mutations --------------------

# Filter expnic edits for syn and NS mutations
exonic_edits <- unique_annotations_Critical %>%
  filter(Func.refGene == "exonic") %>%
  group_by(ExonicFunc.refGene) %>%
  summarize(n_exonic_edits = n())

# Inset plot of Syn/NS
inset_plot <- ggplot(exonic_edits, aes(x = ExonicFunc.refGene, y = n_exonic_edits, fill = ExonicFunc.refGene)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Exonic Function",
    y = "Number of Edits"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),
        legend.position = "none" )
        

# Combine genomic region and inset plot of syn/NS
combined_plot_region_syn_NS_critical <- ggdraw() +
  draw_plot(Unique_critical_genomic_region) +  
  draw_plot(inset_plot, x = 0.6, y = 0.55, width = 0.35, height = 0.35)  
ggsave("Unique_critical_genomic_region_syn_NS.png", width = 10, height = 10)

#========= genomic region_non_critical_unique=======================================

# Group by genomic region and count the number of unique CHROM_POS 
edit_counts_NC <- unique_annotations_Non_critical %>%
  group_by(Func.refGene) %>%
  summarize(n_edits = n_distinct(CHROM_POS)) %>%
  arrange(desc(n_edits)) 

# Convert genomic region to a factor 
edit_counts_NC$Func.refGene <- factor(edit_counts_NC$Func.refGene, levels = edit_counts_NC$Func.refGene[order(edit_counts_NC$n_edits, decreasing = TRUE)])

# Plot
Unique_Non_critical_genomic_region <- ggplot(edit_counts_NC, aes(x = Func.refGene, y = n_edits)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  geom_text(aes(label = n_edits), vjust = -0.3, size = 3) + # Add labels on top of bars
  labs(
    x = "Genomic Region",
    y = "Number of unique edits - critical"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),  
        legend.text = element_text(size = 14),
        legend.position = c(0.9, 0.9),  
        legend.justification = c("right", "top"))

ggsave("Unique_critical_genomic_region.png", width = 8, height = 7)


#---- exonic edits with synonymous and non synonymous mutations --------------------

#filter expnic edits for syn and NS mutations
exonic_edits_NC <- unique_annotations_Non_critical %>%
  filter(Func.refGene == "exonic") %>%
  group_by(ExonicFunc.refGene) %>%
  summarize(n_exonic_edits = n())

# Inset plot of Syn/NS
inset_plot <- ggplot(exonic_edits_NC, aes(x = ExonicFunc.refGene, y = n_exonic_edits, fill = ExonicFunc.refGene)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Exonic Function",
    y = "Number of Edits"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_text(size = 16),
        legend.position = "none" )


# Combine genomic region and inset plot of syn/NS
combined_plot_region_syn_NS_Non_critical <- ggdraw() +
  draw_plot(Unique_Non_critical_genomic_region) +  
  draw_plot(inset_plot, x = 0.6, y = 0.55, width = 0.35, height = 0.35)  
ggsave("Unique__Non_critical_genomic_region_syn_NS.png", width = 10, height = 10)

#=========== over representation analysis =======================================
#https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
#https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/
#https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(DOSE) # cite this - If you use DOSE in published research, please cite:
# Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015, 31(4):608-609
library(enrichplot)

setwd("G:\\Manuscript_2_rerunning_new_SG\\Unique_edits_from_AIDD_VCF_files_annotated_with_Annovar")

# Unique list of genes for both critical and non-critical
critical_genes <- read.csv("unique_sites_annotated_Critical.csv")
gene_list_critical <- unique(critical_genes$Gene.refGene)

noncritical_genes <- read.csv("unique_sites_annotated_Non_critical.csv")
gene_list_noncritical <- unique(noncritical_genes$Gene.refGene)

# Convert to ENTREZ ids
gene_df_critical <- bitr(gene_list_critical, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) # note 14.05% of input gene IDs are fail to map...
gene_df_noncritical <- bitr(gene_list_noncritical, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db) # 13.77% of input gene IDs are fail to map...

# Extract just the ENTREZ IDs
gene_ids_critical <- gene_df_critical$ENTREZID
gene_ids_noncritical <- gene_df_noncritical$ENTREZID

# Combined gene list
gene_lists <- list(
  Critical = gene_ids_critical,
  NonCritical = gene_ids_noncritical
)

# Run ORA - GO-BP
compare_res <- compareCluster(
  geneCluster   = gene_lists,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  keyType       = "ENTREZID",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)
head(compare_res)
?dotplot
both_plots_dot <- dotplot(
  compare_res,
  x = "Cluster",
  color = "p.adjust",
  showCategory = 10,
  split = NULL,
  font.size = 12,
  title = "GO-BP",
  by = "geneRatio",
  size = NULL,
  includeAll = TRUE,
  label_format = 30
)

ggsave("compareCluster_GO_results.png", plot = both_plots_dot, width = 12, height = 10, dpi = 300)

# Save to CSV
write.csv(compare_res, file = "compareCluster_GO_results.csv", row.names = FALSE)

# Run ORA - reactome
compare_res_reactome <- compareCluster(
  geneCluster   = gene_lists,
  fun           = "enrichPathway",
  organism      = "human",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)
# write plot and save the above results
write.csv(compare_res_reactome@compareClusterResult, file = "compare_res_reactome.csv", row.names = FALSE)

reactome_plot <- dotplot(
  compare_res_reactome,
  x = "Cluster",
  color = "p.adjust",
  showCategory = 10,
  split = NULL,
  font.size = 12,
  by = "geneRatio",
  size = NULL,
  includeAll = TRUE,
  label_format = 30
) +
  ggtitle("Reactome Pathway")  

ggsave("compareCluster_Reactome_results.png", plot = reactome_plot, width = 10, height = 10, dpi = 300)

#-----------------------------------------------------------------------------------
## Alu Editing Index
#------------------------------------------------------------------------------------
# Packages
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("G:\\Manuscript_2_rerunning_new_SG\\AEI")
df <- read.csv("EditingIndex.csv")

# Summary stats (mean and SE)
summary_df <- df %>%
  group_by(Disease_severity) %>%
  summarise(
    mean_index = mean(A2GEditingIndex, na.rm = TRUE),
    se = sd(A2GEditingIndex, na.rm = TRUE) / sqrt(n())
  )

# Factor disease severity
df$Disease_severity <- factor(df$Disease_severity, levels = c("Non-critical", "Critical"))
unique(df$Disease_severity)

# Convert Disease_severity to a factor (again, just in case)
df$Disease_severity <- factor(df$Disease_severity, levels = c("Non-critical", "Critical"))

# Reshape from wide to long format:
df_long <- df %>%
  pivot_longer(
    cols = c(A2CEditingIndex, A2GEditingIndex, A2TEditingIndex, C2AEditingIndex, C2GEditingIndex, C2TEditingIndex),
    names_to = "EditingType",
    values_to = "EditingValue"
  )

# Factor EditingType 
df_long <- df %>%
  pivot_longer(
    cols = c(A2CEditingIndex, A2GEditingIndex, A2TEditingIndex, C2AEditingIndex, C2GEditingIndex, C2TEditingIndex),
    names_to = "EditingType",
    values_to = "EditingValue"
  )

df_long$EditingType <- factor(df_long$EditingType,
                              levels = c("A2CEditingIndex", "A2GEditingIndex", "A2TEditingIndex",
                                         "C2AEditingIndex", "C2GEditingIndex", "C2TEditingIndex"),
                              labels = c("A-to-C", "A-to-G", "A-to-T",
                                         "C-to-A", "C-to-G", "C-to-T"))

# Plot 
ggplot(df_long, aes(x = EditingType, y = EditingValue, fill = Disease_severity)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c("Non-critical" = "darkblue", "Critical" = "tomato")) +
  labs(x = "", y = " Alu Editing Index", fill = "Disease Severity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 14))
    legend.title = element_blank(),
    legend.position = c(0.9, 0.9))

ggsave("AEI.png", width = 8, height = 5)

