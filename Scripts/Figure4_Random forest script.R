#------------------------------------------------
# Creating the input matrix for RF analysis
#------------------------------------------------

setwd("G:\\Manuscript_2_rerunning_new_SG\\Random_forest")

# ----------------------------------
# Libraries
# ----------------------------------
library(randomForest)
library(caret)
library(pROC)
library(dplyr)
library(ggplot2)

# 1. List of differentially edited sites
differential_df <- read.csv("G:\\Manuscript_2_rerunning_new_SG\\Random_forest\\edgeR_ATGC_results_FDR0_05_FC0.05_annotated.csv")
differential_df <- differential_df %>%
  mutate(Gene_coordinate = paste0(snpEff_GeneName, "_", CHROM_POS))

# 2. Vector of all the edited sites
gene_coords <- sort(unique(differential_df$Gene_coordinate))

# 3.List of sample files
file_list <- list.files(path = "G:\\Manuscript_2_rerunning_new_SG\\Random_forest\\Jacusa_filtered_4_ADAR_Editing_Level_Added\\", pattern = "\\.csv$", full.names = TRUE)

# 4. Initialize empty matrix with rows = number of files, cols = gene_coords
editing_matrix <- matrix(NA, nrow = length(file_list), ncol = length(gene_coords))
rownames(editing_matrix) <- tools::file_path_sans_ext(basename(file_list))
colnames(editing_matrix) <- gene_coords

# 5. Fill matrix
for (i in seq_along(file_list)) {
  file <- file_list[i]
  df <- read.csv(file)
  
  df <- df %>%
    mutate(Gene_coordinate = paste0(snpEff_GeneName, "_", CHROM_POS)) %>%
    filter(Gene_coordinate %in% gene_coords)
  
  # Map editing levels to corresponding columns
  matched_coords <- df$Gene_coordinate
  matched_values <- df$editing_level
  
  col_indices <- match(matched_coords, gene_coords)
  editing_matrix[i, col_indices] <- matched_values
}

# 6. Convert to data.frame and replace NA with 0 
editing_df <- as.data.frame(editing_matrix)
editing_df[is.na(editing_df)] <- 0  # if desired
# View(editing_df)
condition <- c(rep("non_critical", 23), rep("critical", 46))
editing_df$condition <- condition

write.csv(editing_df, "editing_level_matrix_for_RF.csv", row.names = TRUE)


#----------------------------------------------------------------
# Predicting whether the sample belongs to 'critical' condition
#----------------------------------------------------------------

set.seed(42)

# ----------------------------------
# 1. Data for RF
# ----------------------------------
setwd("G:/Manuscript_2_rerunning_new_SG/Random_forest")
editing_data <- read.csv("editing_level_matrix_for_RF.csv", row.names = 1)
editing_data$condition <- factor(editing_data$condition)
editing_data$condition <- relevel(editing_data$condition, ref = "non_critical")

# ----------------------------------
# 2. Splitting training and test data
# ----------------------------------
set.seed(42)

trainIndex <- createDataPartition(editing_data$condition, p = 0.7, list = FALSE)
trainData <- editing_data[trainIndex, ]
testData  <- editing_data[-trainIndex, ]

# ----------------------------------
# 3. Training the RF model
# ----------------------------------
rf_model <- randomForest(condition ~ ., data = trainData, importance = TRUE,
                         ntree = 500, proximity = TRUE)

# ----------------------------------
# 4. Evaluate the model- Train and Test ROC, AUC
# ----------------------------------
## Training
train_prob <- predict(rf_model, newdata = trainData, type = "prob")
train_roc  <- roc(trainData$condition, train_prob[, "critical"])
train_auc  <- auc(train_roc)

## Testing
test_prob <- predict(rf_model, newdata = testData, type = "prob")
test_roc  <- roc(testData$condition, test_prob[, "critical"])
test_auc  <- auc(test_roc)

## Print AUCs
# print(paste("Training AUC:", round(train_auc, 4)))
# print(paste("Test AUC:", round(test_auc, 4)))

## ROC plots
png("Train_ROC.png")
plot(train_roc, main = "ROC Curve - Training Set", col = "blue", lwd = 2)
dev.off()

png("Test_ROC.png")
plot(test_roc, main = "ROC Curve - Test Set", col = "darkred", lwd = 2)
dev.off()

# ----------------------------------
# 5. Feature importance plot
# ----------------------------------
importance_df <- data.frame(
  Feature = rownames(rf_model$importance),
  MeanDecreaseGini = rf_model$importance[, "MeanDecreaseGini"]
) %>% arrange(desc(MeanDecreaseGini)) %>% slice(1:30)


varImpPlot(rf_model,
           type = 1,               # 1 = accuracy decrease, 2 = Gini decrease
           n.var = 30,            
           main = "Variable Importance (Top 30)")

importance_df <- importance_df %>% arrange(desc(MeanDecreaseGini))

top_features <- importance_df[1:min(30, nrow(importance_df)), ]



p_imp <- ggplot(importance_df, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "steelblue") + coord_flip() 

ggsave("Feature_Importance_Top30.png", plot = p_imp, width = 8, height = 6)

# ----------------------------------
# 6. Per site AUC calculation and plotting
# ----------------------------------
features <- setdiff(names(editing_data), "condition")
site_auc <- sapply(features, function(f) {
  r <- tryCatch(roc(editing_data$condition, editing_data[[f]], levels = levels(editing_data$condition)),
                error = function(e) return(NA))
  if (!is.na(r)[1]) auc(r) else NA
})

auc_df <- data.frame(Feature = names(site_auc), AUC = site_auc)
auc_df <- auc_df[!is.na(auc_df$AUC), ]
top_auc <- auc_df %>% arrange(desc(AUC)) %>% head(11)
write.csv(top_auc, "top_auc.csv", row.names = FALSE)

p_auc <- ggplot(top_auc, aes(x = reorder(Feature, AUC), y = AUC)) +
  geom_col(fill = "darkgreen") + coord_flip() +
  labs(title = "Top 0 Editing Sites by AUC", x = "Editing Site", y = "AUC") +
  theme_minimal()

ggsave("Top_AUC>0.75_AUC_Features.png",  width = 8, height = 6)
ggsave("G:/Manuscript_2_rerunning_new_SG/Random_forest/Top_AUC_gt_0.75_AUC_Features.png", width = 8, height = 6)


