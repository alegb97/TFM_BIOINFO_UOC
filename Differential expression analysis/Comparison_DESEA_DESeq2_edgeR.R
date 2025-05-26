library(VennDiagram)

setwd("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/featureCounts")

# Crear carpeta de salida para los resultados
output_dir <- "GSEA/OBS_Comparison_DESeq2_edgeR"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cargar lista de genes comunes entre modelos en DESEq2 y edgeR 

DESeq2<-read.csv("GSEA/OBS_Comparison_DESeq2/Common_Genes.csv")
edgeR<-read.csv("GSEA/OBS_Comparison_edgeR/Common_Genes.csv")

# Encontrar los genes comunes

genes_comunes <- intersect(DESeq2$x, edgeR$x)

# Guardar los genes comunes 

write.csv(genes_comunes, "genes_comunes.csv", row.names = FALSE)

# Crear diagrama de Venn
tiff("GSEA/OBS_Comparison_DESeq2_edgeR/VennDiagram2.tiff", width = 10, height = 8, units = "in", res = 300)
venn.plot <- draw.pairwise.venn(
  area1 = length(DESeq2$x), area2 = length(edgeR$x), cross.area = length(genes_comunes),
  category = c("DESeq2", "edgeR"),
  fill = c("blue", "red"), alpha = 0.5, cex = 4, cat.cex = 4, cat.pos = c(-20,20)
)
grid.draw(venn.plot)
dev.off()
