# Cargar librerías necesarias
library(DESeq2)
library(ggplot2)
library(pheatmap)

setwd("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/Human_Datasets/Liver/")

# Crear carpeta de salida para los resultados
output_dir <- "GSEA/DESeq2"
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

# Cambiamos nombres columnas
colnames(counts) <- sample_names

# PASO 2: Definir el diseño experimental
colData <- data.frame(
  row.names = sample_names,
  condition = c(rep("Normal-weight", 14), rep("NAFL", 15))
)

# Convertir a factor
colData$condition <- factor(colData$condition, levels = c("Normal-weight", "NAFL"))

# PASO 3: Crear objeto DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)

# Filtrar genes con baja expresión
dds <- dds[rowSums(counts(dds)) > 10, ]

# PASO 4: Ejecutar DESeq2
dds <- DESeq(dds)

# PASO 5: Obtener resultados
res <- results(dds, contrast = c("condition", "Normal-weight", "NAFL"))
res <- res[order(res$padj), ]

# Guardar tabla de resultados en CSV
write.csv(as.data.frame(res), file.path(output_dir, "DESeq2_results_NAFL.csv"))

# Guardar tabla con estadísticos generales
write.csv(as.data.frame(res$stat), file.path(output_dir, "DESeq2_stats_NAFL.csv"))

# Guardar datos de DEGs

DEGs<-subset(as.data.frame(res), (res$log2FoldChange>=1 | res$log2FoldChange<=-1) & res$padj<0.05)
write.csv(DEGs,file.path(output_dir, "DESeq2_results_DEGs_NAFL.csv"))

# PASO 6: Generar gráficos y guardarlos en formato .tiff con 600 dpi

# Histograma de valores ajustados de p
tiff(file.path(output_dir, "histogram_padj_NAFL.tiff"), width = 8, height = 6, units = "in", res = 600)
ggplot(as.data.frame(res), aes(x=padj)) +
  geom_histogram(bins=50, fill="steelblue", color="black") +
  theme_minimal() +
  ggtitle("Distribution of Adjusted p-values") +
  xlab("Adjusted p-value") + 
  ylab("Frequency")
dev.off()

# Volcano plot con colores diferenciados
res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
res$padj[is.na(res$padj)] <- 1

# Definir puntos de corte
logFC_cutoff <- 1  # Punto de corte para log2FoldChange
padj_cutoff <- 0.05  # Punto de corte para p-ajustado

# Crear columna de color
res$color <- "black"
res$color[res$log2FoldChange > logFC_cutoff & res$padj < padj_cutoff] <- "red"
res$color[res$log2FoldChange < -logFC_cutoff & res$padj < padj_cutoff] <- "blue"

# Generar Volcano Plot
tiff(file.path(output_dir, "volcano_plot_NAFL2.tiff"), width = 8, height = 6, units = "in", res = 600)
ggplot(as.data.frame(res), aes(x=log2FoldChange, y=-log10(padj), color=color)) +
  geom_point(alpha=1) +
  scale_color_identity() +
  theme_minimal() +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted p-value") +
  ggtitle("Volcano Plot: Normal-weight vs NAFL") +
  geom_hline(yintercept=-log10(padj_cutoff), linetype="dashed", color="red") +  # Línea horizontal en p-adj 0.05
  geom_vline(xintercept=c(-logFC_cutoff, logFC_cutoff), linetype="dashed", color="black")+ # Líneas verticales para logFC
  theme(axis.title.x =element_text(size = 20), axis.title.y = element_text(size = 20),axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15), plot.title = element_text(size = 20))
dev.off()

# Heatmap de los genes más diferencialmente expresados
topGenes <- rownames(res)[1:50]
heatmap_data <- assay(vst(dds, blind=FALSE))[topGenes,]
tiff(file.path(output_dir, "heatmap_NAFL.tiff"), width = 8, height = 6, units = "in", res = 600)
pheatmap(heatmap_data, cluster_rows=TRUE, cluster_cols=TRUE, show_rownames=FALSE,
         main = "Top 50 Differentially Expressed Genes",
         labels_col = sample_names)
dev.off()

# PASO 7: Análisis de Componentes Principales (PCA)
rld <- rlog(dds, blind=TRUE)  # Transformación logarítmica para PCA
pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Corregir nombres de muestras en el PCA
pcaData$name <- as.character(pcaData$name) # Asegura que sean caracteres y no factores
pcaData$condition <- factor(pcaData$condition, levels = c("Normal-weight", "NAFL"))


tiff(file.path(output_dir, "PCA_plot_NAFL2.tiff"), width = 8, height = 6, units = "in", res = 600)
ggplot(pcaData, aes(x=PC1, y=PC2, color=condition, label=name)) +
  geom_point(size=8) +
  theme_minimal() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("Principal Component Analysis (PCA)") +
  geom_text(vjust=-1.5, hjust=0.5, size=3)+ # Ajusta la posición de las etiquetas
  theme(axis.title.x =element_text(size = 20), axis.title.y = element_text(size = 20),axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 15), plot.title = element_text(size = 20),legend.title =element_text(size = 15) ,legend.text = element_text(size = 15))

dev.off()

