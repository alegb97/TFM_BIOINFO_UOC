library(RRHO2)
library(readr)

setwd("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/featureCounts")

# Definir nombres de salida
output_dir <- "GSEA/OBS_Comparison_RRHO2/edgeR"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cargar listas de genes diferencialmente expresados
res_dio <- read.csv("GSEA/DIO/edgeR/edgeR_results.csv", row.names = 1)
res_zucker <- read.csv("GSEA/Zucker/edgeR/edgeR_results.csv", row.names = 1)

# Crear listas ordenadas para RRHO2
gene_dio <- data.frame(
  Gene = rownames(res_dio), 
  Score = -log10(res_dio$FDR) * sign(res_dio$logFC)
)

gene_zucker <- data.frame(
  Gene = rownames(res_zucker), 
  Score = -log10(res_zucker$FDR) * sign(res_zucker$logFC)
)

# Filtrar valores NA
gene_dio <- na.omit(gene_dio)
gene_zucker <- na.omit(gene_zucker)

# Filtrar los genes comunes en los DataFrames
common_genes <- intersect(gene_dio$Gene, gene_zucker$Gene)
gene_dio <- gene_dio[gene_dio$Gene %in% common_genes, ]
gene_zucker <- gene_zucker[gene_zucker$Gene %in% common_genes, ]

# Ordenar las listas de mayor a menor puntuación
gene_dio <- gene_dio[order(-gene_dio$Score), ]
gene_zucker <- gene_zucker[order(-gene_zucker$Score), ]

# Verificar que tienen los mismos genes
identical(gene_dio$Gene, gene_zucker$Gene)

# Ejecutar RRHO2
rrho_result <- RRHO2_initialize(gene_dio, gene_zucker, labels = c("OBS-DIO", "OBS-GEN"), log10.ind = TRUE)

# Guardar heatmap
tiff(file.path(output_dir, "RRHO2_Heatmap.tiff"), width = 8, height = 6, units = "in", res = 600)
RRHO2_heatmap(rrho_result)
dev.off()

# Guardar diagrama de Venn Down-regulated-Down-regulated
tiff(file.path(output_dir, "RRHO2_VennDiagram_DD.tiff"), width = 8, height = 6, units = "in", res = 600)
RRHO2_vennDiagram(rrho_result, type = "dd")
dev.off()

# Guardar diagrama de Venn Up-regulated-Up-regulated
tiff(file.path(output_dir, "RRHO2_VennDiagram_UU.tiff"), width = 8, height = 6, units = "in", res = 600)
RRHO2_vennDiagram(rrho_result, type = "uu")
dev.off()

# Guardar diagrama de Venn Up-regulated-Down-regulated
tiff(file.path(output_dir, "RRHO2_VennDiagram_UD.tiff"), width = 8, height = 6, units = "in", res = 600)
RRHO2_vennDiagram(rrho_result, type = "ud")
dev.off()

# Guardar diagrama de Venn Down-regulated-Up-regulated
tiff(file.path(output_dir, "RRHO2_VennDiagram_DU.tiff"), width = 8, height = 6, units = "in", res = 600)
RRHO2_vennDiagram(rrho_result, type = "du")
dev.off()

# Guardar genes sobreexpresados e infraexpresados en ambos modelos

write_csv(as.data.frame(rrho_result$genelist_uu$gene_list_overlap_uu),"GSEA/OBS_Comparison_RRHO2/edgeR/upregulated_DIO_Zucker.csv")
write_csv(as.data.frame(rrho_result$genelist_dd$gene_list_overlap_dd),"GSEA/OBS_Comparison_RRHO2/edgeR/downregulated_DIO_Zucker.csv")



