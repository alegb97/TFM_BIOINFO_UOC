# Cargar librerías necesarias
library(edgeR)
library(ggplot2)
library(pheatmap)

# Definir directorio de trabajo
setwd("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/featureCounts")

# Crear carpeta de salida
output_dir <- "GSEA/DIO/edgeR"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cargar datos de featureCounts
counts <- read.csv("featureCounts_DIO_results.csv", row.names = 1)

# Definir nombres de las muestras (13 semanas)
samples_13w <- c("DIO_S1_sorted.bam", "DIO_S2_sorted.bam", "DIO_S3_sorted.bam",
                 "DIO_S4_sorted.bam", "DIO_S5_sorted.bam", "DIO_S6_sorted.bam",
                 "DIO_S7_sorted.bam", "DIO_S8_sorted.bam", "DIO_S9_sorted.bam",
                 "DIO_S10_sorted.bam", "DIO_S11_sorted.bam", "DIO_S12_sorted.bam")

# Seleccionar solo esas columnas
counts <- counts[, samples_13w]

# Limpiar nombres de las muestras
sample_names <- gsub("_sorted.bam", "", samples_13w)
colnames(counts) <- sample_names

# Definir diseño experimental
group <- factor(c(rep("CT-DIO", 6), rep("OBS-DIO", 6)))

# Crear objeto DGEList
dge <- DGEList(counts = counts, group = group)

# Filtrar genes con muy baja expresión
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
write.csv(res_table, file.path(output_dir, "edgeR_results.csv"))

# Guardar datos de DEGs

DEGs<-subset(res_table, (res_table$logFC>=1 | res_table$logFC<=-1) & res_table$FDR<0.05)
write.csv(DEGs,file.path(output_dir, "edgeR_results_DEGs.csv") )

# **Volcano Plot**
logFC_cutoff <- 1
fdr_cutoff <- 0.05

res_table$color <- "black"
res_table$color[res_table$logFC > logFC_cutoff & res_table$FDR < fdr_cutoff] <- "red"
res_table$color[res_table$logFC < -logFC_cutoff & res_table$FDR < fdr_cutoff] <- "blue"

tiff(file.path(output_dir, "volcano_plot.tiff"), width = 8, height = 6, units = "in", res = 600)
ggplot(res_table, aes(x = logFC, y = -log10(FDR), color = color)) +
  geom_point(alpha = 1) +
  scale_color_identity() +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 FDR") +
  ggtitle("Volcano Plot: CT-DIO vs OBS-DIO") +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black")+
  theme(axis.title.x =element_text(size = 20), axis.title.y = element_text(size = 20),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15), plot.title = element_text(size = 20))

dev.off()

# **Heatmap de los 50 genes más significativos**
topGenes <- rownames(res_table)[1:50]
heatmap_data <- cpm(dge, log = TRUE)[topGenes, , drop = FALSE]

tiff(file.path(output_dir, "heatmap.tiff"), width = 8, height = 6, units = "in", res = 600)
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

tiff(file.path(output_dir, "PCA_plot2.tiff"), width = 8, height = 6, units = "in", res = 600)
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 8) +
  theme_minimal() +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA: CT-DIO vs OBS-DIO") +
  geom_text(vjust = -1.5, hjust = 0.5, size = 3)+
  theme(axis.title.x =element_text(size = 20), axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 15), plot.title = element_text(size = 20),legend.title =element_text(size = 15) ,legend.text = element_text(size = 15))
dev.off()
