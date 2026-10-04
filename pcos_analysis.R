
dir.create("C:/bioinformatics_projects/pcos_shap", recursive = TRUE, 
           showWarnings = FALSE)
setwd ("C:/bioinformatics_projects/pcos_shap")

dir.create("data/GSE138518", recursive = TRUE, showWarnings = FALSE)

library(GEOquery)
gse_files <- getGEOSuppFiles("GSE138518", makeDirectory = FALSE,
                             baseDir = "data/GSE138518")
gse_files

# Step 2: Install readxl and inspect the file

install.packages("readxl")
library(readxl)

sheets <- excel_sheets("data/GSE138518/GSE138518_RNA.xlsx")
sheets

preview <- read_excel("data/GSE138518/GSE138518_RNA.xlsx", sheet = 1,
                       n_max = 5)
dim(preview)
preview[ , 1:5]

# Step 3: See all column names to identify which samples are PCOS vs Control

colnames(preview)

# Step 4: Load the full dataset and build samole metadata

full_data <- read_excel("data/GSE138518/GSE138518_RNA.xlsx", sheet = 1)
dim(full_data)

raw_counts <- full_data[, c("ENSEMBL", "SYMBOL", "N20", "N21", "N25", 
                            "P14", "P15", "P16")]
raw_counts

# Step 5: Build the DESeq2 object and normalize

count_matrix <- as.data.frame(raw_counts[, c("N20", "N21", "N25", 
                                             "P14", "P15", "P16")])
rownames(count_matrix) <- raw_counts$ENSEMBL
count_matrix <- round(count_matrix)

gene_lookup <- data.frame(ENSEMBL = raw_counts$ENSEMBL, SYMBOL = raw_counts$SYMBOL)
head(gene_lookup)

sample_info <- data.frame(
  condition = c("control", "control", "control", "PCOS", "PCOS", "PCOS"),
  row.names = colnames(count_matrix)
)
sample_info

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

dds_DESeq2 <- DESeq(dds)

# Step 6: Export normalized expression for the classifier

norm_counts <- counts(dds_DESeq2, normalized = TRUE)
log_counts <- log2(norm_counts + 1)

export_data <- as.data.frame(t(log_counts))
export_data$condition <- sample_info$condition

write.csv(export_data, "pcos_classifier_input.csv", row.names = TRUE)

# Step 7: Also export the gene lookup table

write.csv(gene_lookup, "gene_lookup.csv", row.names = FALSE)



# Step 18: Map the SHAP selected genes to STRING

library(STRINGdb)

string_db <- STRINGdb$new(version = "12.0", species = 9606, 
                          score_threshold = 400, input_directory = "")

shap_genes <- data.frame(gene =c("FPR2", "VNN3", "E2F8", "CCL25",
                                 "USH2A", "P3H3", "AVPR1B", "SLC35F1"))

mapped_shap_genes <- string_db$map(shap_genes, "gene", 
                                   removeUnmappedRows = TRUE)
mapped_shap_genes

# Step 19: Check how the 8 genes interract with each other

string_db$plot_network(mapped_shap_genes$STRING_id)
