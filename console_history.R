getwd()
list.files()
list.files("data")
pkgs <- c("tidyverse","DESeq2","pheatmap","glmnet","survival",
"clusterProfiler","org.Hs.eg.db","ggplot2")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install) > 0) install.packages(to_install)
if(!"DESeq2" %in% rownames(installed.packages())) {
if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
}
rna <- read.delim("data/data_mrna_seq_v2_rsem.txt",
check.names = FALSE, stringsAsFactors = FALSE)
cna <- read.delim("data/data_cna.txt",
check.names = FALSE, stringsAsFactors = FALSE)
clin <- read.delim("data/data_clinical_patient.txt",
check.names = FALSE, stringsAsFactors = FALSE)
dim(rna); dim(cna); dim(clin)
head(colnames(rna), 10)
head(colnames(cna), 10)
head(colnames(clin), 10)
head(rna[,1:3])
head(cna[,1:3])
head(clin[,1:3])
clin_clean <- clin %>%
filter(`#Patient Identifier` != "#Patient Identifier") %>%
filter(`#Patient Identifier` != "PATIENT_ID")
library(tidyverse)
clin_clean <- clin %>%
filter(`#Patient Identifier` != "#Patient Identifier") %>%
filter(`#Patient Identifier` != "PATIENT_ID")
head(clin_clean[, 1:3])
dim(clin_clean)
clin_clean <- clin_clean[-(1:3), ]
head(clin_clean[, 1:3])
dim(clin_clean)
rna_patients <- substr(colnames(rna)[-(1:2)], 1, 12)
cna_patients <- substr(colnames(cna)[-(1:2)], 1, 12)
clin_patients <- clin_clean$`#Patient Identifier`
common_patients <- Reduce(intersect, list(rna_patients, cna_patients, clin_patients))
length(common_patients)
# Subset RNA
rna_keep <- rna_patients %in% common_patients
rna_mat <- rna[, c(TRUE, TRUE, rna_keep)]
colnames(rna_mat)[-(1:2)] <- substr(colnames(rna_mat)[-(1:2)], 1, 12)
# Subset CNA
cna_keep <- cna_patients %in% common_patients
cna_mat <- cna[, c(TRUE, TRUE, cna_keep)]
colnames(cna_mat)[-(1:2)] <- substr(colnames(cna_mat)[-(1:2)], 1, 12)
all(colnames(rna_mat)[-(1:2)] == colnames(cna_mat)[-(1:2)])
erbb2_cna <- cna_mat[cna_mat$Hugo_Symbol == "ERBB2", ]
dim(erbb2_cna)  # should be 1 row
metadata <- data.frame(
patient_id = colnames(erbb2_cna)[-(1:2)],
ERBB2_cna = as.numeric(erbb2_cna[1, -(1:2)])
)
metadata$ERBB2_amp <- metadata$ERBB2_cna > 0
table(metadata$ERBB2_amp)
expr_mat <- rna_mat[, -(1:2)]
rownames(expr_mat) <- rna_mat$Hugo_Symbol
expr_mat <- rna_mat[, -(1:2)]
gene <- rna_mat$Hugo_Symbol
# Replace blank / NA gene symbols with Entrez IDs
gene[is.na(gene) | gene == ""] <- paste0("ENTREZ_", rna_mat$Entrez_Gene_Id[is.na(gene) | gene == ""])
# Make duplicated gene symbols unique by appending .1 .2 etc
gene <- make.unique(gene)
rownames(expr_mat) <- gene
any(duplicated(rownames(expr_mat)))
sum(rownames(expr_mat) == "" | is.na(rownames(expr_mat)))
dim(expr_mat)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(metadata, "output/tables/metadata_matched.csv", row.names = FALSE)
write.csv(expr_mat, "output/tables/expression_matched.csv")
dir.create("Output_tables", showWarnings = FALSE)
dir.create("Output_figures", showWarnings = FALSE)
write.csv(metadata, "Output_tables/metadata_matched.csv", row.names = FALSE)
write.csv(expr_mat, "Output_tables/expression_matched.csv")
library(tidyverse)
# Ensure output folders exist (your naming style)
dir.create("Output_tables", showWarnings = FALSE)
dir.create("Output_figures", showWarnings = FALSE)
# --- Rebuild expression matrix safely (unique rownames) ---
expr_mat <- rna_mat[, -(1:2)]
gene <- rna_mat$Hugo_Symbol
# Replace blank/NA symbols with Entrez IDs
bad <- is.na(gene) | gene == ""
gene[bad] <- paste0("ENTREZ_", rna_mat$Entrez_Gene_Id[bad])
# Make duplicates unique (adds .1, .2 etc)
gene <- make.unique(gene)
rownames(expr_mat) <- gene
# --- Save checkpoint files ---
write.csv(metadata, "Output_tables/metadata_matched.csv", row.names = FALSE)
write.csv(expr_mat, "Output_tables/expression_matched.csv")
# --- Check ERBB2 amplification counts ---
print(table(metadata$ERBB2_amp))
# --- Check whether expression values are integers (DESeq2 needs integers) ---
x <- as.numeric(expr_mat[1:1000, 1])
x <- x[is.finite(x)]
frac_decimals <- mean(x %% 1 != 0)
print(frac_decimals)
summary(x)
# Round to integers for DESeq2 (required)
counts_mat <- round(as.matrix(expr_mat))
counts_mat[counts_mat < 0] <- 0
# Align metadata row order to counts columns
meta <- metadata %>%
mutate(ERBB2_amp = factor(ERBB2_amp, levels = c(FALSE, TRUE))) %>%
column_to_rownames("patient_id")
counts_mat <- counts_mat[, rownames(meta)]
# Prefilter low-expression genes (speeds up + improves stability)
keep <- rowSums(counts_mat >= 10) >= 10
counts_filt <- counts_mat[keep, ]
dim(counts_filt)
write.csv(counts_filt, "Output_tables/counts_for_deseq2_filtered.csv")
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
countData = counts_filt,
colData   = meta,
design    = ~ ERBB2_amp
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("ERBB2_amp", "TRUE", "FALSE"))
res_tbl <- as.data.frame(res) %>%
rownames_to_column("gene") %>%
arrange(padj)
# Save full DE table
write.csv(res_tbl, "Output_tables/DESeq2_results_all.csv", row.names = FALSE)
# Significant DE genes (FDR < 0.05)
sig <- res_tbl %>% filter(!is.na(padj), padj < 0.05)
write.csv(sig, "Output_tables/DESeq2_results_significant_FDR05.csv", row.names = FALSE)
# Top 10 by absolute log2FC (among significant)
top10 <- sig %>%
arrange(desc(abs(log2FoldChange))) %>%
slice(1:10)
# Convert DESeq2 results to plain data.frame (no Rle weirdness)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
# Ensure numeric types (extra safety)
res_df$log2FoldChange <- as.numeric(res_df$log2FoldChange)
res_df$padj <- as.numeric(res_df$padj)
# Order by adjusted p-value
res_df <- res_df[order(res_df$padj), ]
# Save full DE table
write.csv(res_df, "Output_tables/DESeq2_results_all.csv", row.names = FALSE)
# Significant genes (FDR < 0.05)
sig <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ]
write.csv(sig, "Output_tables/DESeq2_results_significant_FDR05.csv", row.names = FALSE)
# Top 10 by absolute log2FC among significant
sig2 <- sig[order(-abs(sig$log2FoldChange)), ]
top10 <- head(sig2, 10)
write.csv(top10, "Output_tables/Top10_DE_genes_by_abs_log2FC.csv", row.names = FALSE)
n_sig <- nrow(sig)
print(n_sig)
print(top10[, c("gene","log2FoldChange","pvalue","padj")])
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)
dim(vst_mat)
pca <- prcomp(t(vst_mat), scale. = FALSE)
pca_df <- data.frame(
PC1 = pca$x[,1],
PC2 = pca$x[,2],
ERBB2_amp = meta$ERBB2_amp
)
library(ggplot2)
p <- ggplot(pca_df, aes(PC1, PC2, color = ERBB2_amp)) +
geom_point(alpha = 0.7, size = 2) +
labs(
title = "PCA of TCGA Breast Cancer (VST expression)",
x = "PC1",
y = "PC2",
color = "ERBB2 amplified"
) +
theme_minimal()
print(p)
ggsave("Output_figures/PCA_ERBB2_VST.png",
plot = p, width = 7, height = 5)
top50_genes <- res_df$gene[1:50]
length(top50_genes)
heat_mat <- vst_mat[top50_genes, ]
heat_mat_scaled <- t(scale(t(heat_mat)))
annotation_col <- data.frame(
ERBB2_amp = meta$ERBB2_amp
)
rownames(annotation_col) <- rownames(meta)
library(pheatmap)
pheatmap(
heat_mat_scaled,
annotation_col = annotation_col,
show_colnames = FALSE,
show_rownames = TRUE,
fontsize_row = 6,
filename = "Output_figures/Heatmap_Top50_DE_genes.png",
width = 8,
height = 10
)
# Look for likely OS time + status columns
cn <- colnames(clin_clean)
time_candidates <- cn[grepl("overall|survival|os|days|months|time", cn, ignore.case = TRUE)]
event_candidates <- cn[grepl("vital|status|death|deceased|event", cn, ignore.case = TRUE)]
cat("TIME candidates:\n")
print(time_candidates)
cat("\nEVENT/STATUS candidates:\n")
print(event_candidates)
library(tidyverse)
library(survival)
# Build survival dataframe from clinical
surv_df <- clin_clean %>%
transmute(
patient_id = `#Patient Identifier`,
os_months = as.numeric(`Overall Survival (Months)`),
os_status_raw = `Overall Survival Status`
)
# Clean event coding:
# cBioPortal often uses strings like "LIVING" / "DECEASED"
surv_df$event <- ifelse(grepl("DECEASED|DEAD|1", surv_df$os_status_raw, ignore.case = TRUE), 1, 0)
# Keep only patients we have in metadata (and drop missing time)
surv_df <- surv_df %>%
filter(patient_id %in% metadata$patient_id) %>%
filter(!is.na(os_months) & os_months > 0)
# Align to meta order
surv_df <- surv_df %>%
slice(match(rownames(meta), patient_id))
# Build survival dataframe (base-safe)
surv_df <- data.frame(
patient_id = clin_clean$`#Patient Identifier`,
os_months = suppressWarnings(as.numeric(clin_clean$`Overall Survival (Months)`)),
os_status_raw = clin_clean$`Overall Survival Status`,
stringsAsFactors = FALSE
)
# Event coding: DECEASED = 1, otherwise 0 (usually LIVING)
surv_df$event <- ifelse(grepl("DECEASED|DEAD", surv_df$os_status_raw, ignore.case = TRUE), 1, 0)
# Keep only patients in meta and with valid time
surv_df <- surv_df[surv_df$patient_id %in% rownames(meta), ]
surv_df <- surv_df[!is.na(surv_df$os_months) & surv_df$os_months > 0, ]
# Align to meta row order using match (no slice)
idx <- match(rownames(meta), surv_df$patient_id)
surv_df <- surv_df[idx, ]
# Drop any patients missing survival after matching
keep_surv <- !is.na(surv_df$patient_id)
surv_df <- surv_df[keep_surv, ]
# Also align meta accordingly (critical!)
meta_surv <- meta[surv_df$patient_id, , drop = FALSE]
# Sanity checks
stopifnot(all(rownames(meta_surv) == surv_df$patient_id))
summary(surv_df$os_months)
table(surv_df$event)
library(glmnet)
library(survival)
# Pick DE genes: significant then top 500 by padj
de_genes <- sig$gene
de_genes <- de_genes[de_genes %in% rownames(vst_mat)]
de_genes <- head(de_genes, 500)
# X matrix: patients x genes (only patients with survival)
X <- t(vst_mat[de_genes, colnames(vst_mat) %in% surv_df$patient_id, drop = FALSE])
X <- X[surv_df$patient_id, , drop = FALSE]
# Survival object
y <- Surv(time = surv_df$os_months, event = surv_df$event)
set.seed(1)
cvfit <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)
fit <- glmnet(X, y, family = "cox", alpha = 1, lambda = cvfit$lambda.min)
beta <- as.matrix(coef(fit))
sel <- beta[beta[,1] != 0, , drop = FALSE]
sel_df <- data.frame(
gene = rownames(sel),
coef = sel[,1],
row.names = NULL
)
sel_df <- sel_df[order(-abs(sel_df$coef)), ]
nrow(sel_df)
head(sel_df, 20)
write.csv(sel_df, "Output_tables/CoxLASSO_selected_genes.csv", row.names = FALSE)
# Risk score
risk <- as.numeric(predict(fit, newx = X, type = "link"))
surv_df$risk <- risk
surv_df$risk_group <- ifelse(surv_df$risk >= median(surv_df$risk, na.rm = TRUE), "High", "Low")
km_fit <- survfit(Surv(os_months, event) ~ risk_group, data = surv_df)
png("Output_figures/KM_risk_groups.png", width = 900, height = 700)
plot(km_fit, col = c("blue","red"), lwd = 2,
xlab = "Overall Survival (Months)", ylab = "Survival Probability",
main = "Kaplan–Meier: Cox LASSO Risk Groups")
legend("bottomleft", legend = c("High","Low"), col = c("blue","red"), lwd = 2, bty = "n")
dev.off()
write.csv(surv_df, "Output_tables/Survival_table_used.csv", row.names = FALSE)
savehistory("console_history.R")
