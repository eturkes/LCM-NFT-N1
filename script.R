#    This file is part of LCM-NFT-N1.
#    Copyright (C) 2023  Emir Turkes, Martha Foiani, UK DRI at UCL
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Emir Turkes can be contacted at emir.turkes@eturkes.com

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
hub_gene_data_big <- read.csv(file.path("data", "hub_gene_data_big.tsv"), sep = "")
probes_86 <- read.csv(file.path("data", "NFT-probes-86.txt"), sep = "")
probe_names <- hub_gene_data[hub_gene_data$ensembl_gene_id %in% probes_86$initial_alias, ]
probe_names <- unique(probe_names$external_gene_name)

LCM_names <- rownames(results)
keep <- which(LCM_DEPs_raw_data$Protein.Names %in% LCM_names)
LCM_names <- LCM_DEPs_raw_data$Genes[keep]

intersection <- probe_names[probe_names %in% LCM_names]
intersection_all <- intersect(hub_gene_data$external_gene_name, LCM_names)
intersection_all_big <- intersect(hub_gene_data_big$external_gene_name, LCM_names)
intersection_all_big_DE <- hub_gene_data_big$external_gene_name[hub_gene_data_big$DE == "yes"]
intersection_all_big_DE <- intersect(intersection_all_big_DE, LCM_names)
intersection_all_big_DE_data <- hub_gene_data_big[hub_gene_data_big$external_gene_name %in% intersection_all_big_DE, ]

prot_names <- which(LCM_DEPs_raw_data$Genes %in% intersection)
prot_names <- LCM_DEPs_raw_data$Protein.Names[prot_names]

heatmap_mat <- mat[rownames(mat) %in% prot_names, ]
heatmap_mat <- t(apply(heatmap_mat, 1, function (x) ((2 * (x - min(x)) / (max(x) - min(x))) - 1)))
Heatmap(heatmap_mat, cluster_columns = FALSE, cluster_rows = FALSE)
