# Cargar librerías necesarias
library(edgeR)
library(ggplot2)
library(pheatmap)

setwd("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/Human_Datasets/Liver/")

# Crear carpeta de salida para los resultados
output_dir <- "GSEA/edgeR"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# PASO 1: Cargar datos del conteo (hecho con LIMMA) del GSE126848 (Liver dataset)

counts <- read.table("GSE126848_Gene_counts_raw.txt", row.names = 1, header=TRUE, check.names = FALSE)

# Filtramos las muestras de personas con peso normal y obesas
samples<-c("0869","0873","0875","0877","0879","0881","0883","0885","0887","0889","0891","0893","0897","0910",
           "2683","2685","2687","2689","2691","2693","2697","2701","2703","3993","3995","3998","4002","4007","4010")

counts <- counts[, samples]


# Nombres muestras
sample_names<- c("Normal-weight_1","Normal-weight_2","Normal-weight_3","Normal-weight_4","Normal-weight_5","Normal-weight_6","Normal-weight_7","Normal-weight_8","Normal-weight_9",
                 "Normal-weight_10","Normal-weight_11","Normal-weight_12","Normal-weight_13","Normal-weight_14","NAFL_1","NAFL_2","NAFL_3","NAFL_4","NAFL_5","NAFL_6","NAFL_7",
                 "NAFL_8","NAFL_9","NAFL_10","NAFL_11","NAFL_12","NAFL_13","NAFL_14","NAFL_15")

# Cambiamos nombres datasets
colnames(counts) <- sample_names

# Definir diseño experimental
group <- factor(c(rep("Normal-weight", 14), rep("NAFL", 15)))

# Crear objeto DGEList
dge <- DGEList(counts = counts, group = group)

# Filtrar genes con baja expresión
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalizar datos con TMM
dge <- calcNormFactors(dge)

# Estimar dispersión
dge <- estimateDisp(dge)

# Ajuste del modelo y prueba exacta de edgeR
fit <- glmQLFit(dge, design = model.matrix(~group))
res <- glmQLFTest(fit, coef = 2)

# Convertir en un dataframe ordenado
res_table <- topTags(res, n = nrow(res))$table
res_table <- res_table[order(res_table$FDR), ]

# Guardar tabla de resultados
write.csv(res_table, file.path(output_dir, "edgeR_results_NAFL.csv"))

# Guardar datos de DEGs

DEGs<-subset(res_table, (res_table$logFC>=1 | res_table$logFC<=-1) & res_table$FDR<0.05)
write.csv(DEGs,file.path(output_dir, "edgeR_results_DEGs_NAFL.csv") )

# **Volcano Plot**
logFC_cutoff <- 1
fdr_cutoff <- 0.05

res_table$color <- "black"
res_table$color[res_table$logFC > logFC_cutoff & res_table$FDR < fdr_cutoff] <- "red"
res_table$color[res_table$logFC < -logFC_cutoff & res_table$FDR < fdr_cutoff] <- "blue"

tiff(file.path(output_dir, "volcano_plot_NAFL.tiff"), width = 8, height = 6, units = "in", res = 600)
ggplot(res_table, aes(x = logFC, y = -log10(FDR), color = color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 FDR") +
  ggtitle("Volcano Plot: Normal-weight vs NAFL") +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black")
dev.off()

# **Heatmap de los 50 genes más significativos**
topGenes <- rownames(res_table)[1:50]
heatmap_data <- cpm(dge, log = TRUE)[topGenes, , drop = FALSE]

tiff(file.path(output_dir, "heatmap_NAFL.tiff"), width = 8, height = 6, units = "in", res = 600)
pheatmap(heatmap_data, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = FALSE,
         main = "Top 50 Differentially Expressed Genes",
         labels_col = sample_names)
dev.off()

# **Análisis PCA**
pca <- prcomp(t(cpm(dge, log = TRUE)))

pcaData <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  sample = sample_names,
  condition = group
)

tiff(file.path(output_dir, "PCA_plot_NAFL.tiff"), width = 8, height = 6, units = "in", res = 600)
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4) +
  theme_minimal() +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA: Normal-weight vs NAFL") +
  geom_text(vjust = -1.5, hjust = 0.5, size = 3)
dev.off()
