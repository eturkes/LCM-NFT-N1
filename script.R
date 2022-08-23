library(readxl)
library(limma)
library(ComplexHeatmap)

LCM_DEPs_raw_data <- read_excel(file.path("data", "LCM_DEPs raw data.xlsx"))

colnames(LCM_DEPs_raw_data) <- c(
  colnames(LCM_DEPs_raw_data)[1:13], paste(rep("AT8_pos", 3), 1:3, sep = "_"), paste(rep("AT8_neg", 4), 1:4, sep = "_")
)

mat <- as.matrix(LCM_DEPs_raw_data[ , 14:ncol(LCM_DEPs_raw_data)])
rownames(mat) <- LCM_DEPs_raw_data$Protein.Names

design <- model.matrix(~ 0 + sub("_[^_]+$", "", colnames(mat)))
colnames(design) <- unique(sub("_[^_]+$", "", colnames(mat)))
fit <- lmFit(mat, design)
contrast_mat <- makeContrasts(AT8_pos-AT8_neg, levels = design)
cont_fit <- eBayes(contrasts.fit(fit, contrast_mat))
results <- topTable(cont_fit, number = Inf, adjust.method = "none", p.value = 0.05)

heatmap_mat <- mat[rownames(mat) %in% rownames(results)[1:25], ]
heatmap_mat <- t(apply(heatmap_mat, 1, function (x) ((2 * (x - min(x)) / (max(x) - min(x))) - 1)))
Heatmap(heatmap_mat, cluster_columns = FALSE)

hub_gene_data <- read.csv(file.path("data", "hub_gene_data_86_169.tsv"), sep = "")
probes_86 <- read.csv(file.path("data", "NFT-probes-86.txt"), sep = "")
probe_names <- hub_gene_data[hub_gene_data$ensembl_gene_id %in% probes_86$initial_alias, ]
probe_names <- unique(probe_names$external_gene_name)

LCM_names <- sub("_[^_]+$", "", rownames(results))

intersection <- probe_names[probe_names %in% LCM_names]
intersection_all <- unique(
  hub_gene_data$external_gene_name[hub_gene_data$external_gene_name %in% LCM_names]
)
