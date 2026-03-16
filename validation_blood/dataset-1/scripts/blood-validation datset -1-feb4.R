# ======================================================
# BLOOD DATASET 1 VALIDATION
# ======================================================

rm(list=ls())
set.seed(123)

library(pROC)
library(ggpubr)
library(randomForest)
library(caret)

# Final biomarkers
final_genes <- c("AHNAK2","TRIM29","SLC6A14","ITGB4")

# Load data
expr <- read.csv("~/work/r_project_ishu/validation_blood/normalized_expression_matrix_blood.csv", row.names=1, check.names=FALSE)
meta <- read.csv("~/work/r_project_ishu/validation_blood/metadata_blood_clean.csv")
head(meta)
# Convert disease_state to group
meta$group <- meta$disease_state

# Make it a factor
meta$group <- factor(meta$group, levels=c("control","PDAC"))
table(meta$group)

# Extract biomarkers
common_genes <- intersect(final_genes, rownames(expr))
expr_biomarker <- t(expr[common_genes, ])
df <- cbind(meta, expr_biomarker)

cat("Genes used:", common_genes, "\n")

# ---------------- Boxplots ----------------
for(g in common_genes){
  p <- ggboxplot(df, x="group", y=g, color="group", add="jitter") +
    stat_compare_means(method="wilcox.test") +
    ggtitle(paste("Dataset1 Blood Expression:", g))
  print(p)
}

# ---------------- Single Gene ROC ----------------
roc_table <- data.frame(Gene=character(), AUC=numeric())

for(g in common_genes){
  roc_obj <- roc(df$group, df[,g])
  plot(roc_obj, main=paste("ROC Dataset1:", g))
  roc_table <- rbind(roc_table, data.frame(Gene=g, AUC=auc(roc_obj)))
}
write.csv(roc_table, "~/work/r_project_ishu/validation_blood/Dataset1_SingleGene_AUC_feb18.csv", row.names=FALSE)

# ---------------- Multi-Gene ROC ----------------
glm_model <- glm(group ~ ., data=df[,c("group", common_genes)], family="binomial")
pred <- predict(glm_model, type="response")
roc_panel <- roc(df$group, pred)
plot(roc_panel, main="Dataset1 Multi-Gene ROC")
auc(roc_panel)

# ---------------- Random Forest Validation ----------------
trainIndex <- createDataPartition(df$group, p=0.7, list=FALSE)
train <- df[trainIndex,]
test  <- df[-trainIndex,]

rf_model <- randomForest(group ~ ., data=train[,c("group", common_genes)], ntree=500)
pred_prob <- predict(rf_model, test, type="prob")[,2]

roc_rf <- roc(test$group, pred_prob)
plot(roc_rf, main="Dataset1 RF ROC")
auc(roc_rf)

# Confusion matrix
pred_class <- ifelse(pred_prob > 0.5, levels(df$group)[2], levels(df$group)[1])
confusionMatrix(as.factor(pred_class), test$group)
