# Cargar paquetes necesarios
library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Rn.eg.db)
library(biomaRt)
library(enrichplot)
library(ggplot2)

setwd("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/")

# Crear carpetas para resultados
dir.create("Homólogos", showWarnings = FALSE)
dir.create("Homólogos/GO", showWarnings = FALSE)
dir.create("Homólogos/KEGG", showWarnings = FALSE)

# Función de análisis para un modelo de obesidad
analisis_funcional <- function(ruta_csv, modelo) {
  mensaje <- paste("Procesando", modelo, "...")
  cat(mensaje, "\n")
  
  # Leer datos
  data <- read_csv(ruta_csv)
  colnames(data)[1] <- "gene"  # la columna de símbolos no tiene nombre
  
  # Filtrar genes significativos
  data_filtrada <- data %>%
    filter(abs(log2FoldChange) > 1, padj < 0.05) %>%
    filter(!is.na(padj))
  
  # Convertir símbolos a mayúsculas
  genes <- toupper(unique(data_filtrada$gene))
  
  # Conectar con Ensembl para Rattus norvegicus
  rat_mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  
  # Mapear SYMBOL -> ENTREZ
  gene_map <- getBM(
    filters = "external_gene_name",
    attributes = c("external_gene_name", "entrezgene_id"),
    values = genes,
    mart = rat_mart
  )
  
  # Verificar mapeo
  if (nrow(gene_map) == 0) {
    cat("No se pudieron mapear genes para", modelo, "\n")
    return(NULL)
  }
  
  entrez_ids <- unique(na.omit(gene_map$entrezgene_id))
  
  ## Análisis GO
  ego <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ## Análisis KEGG
  ekegg <- enrichKEGG(
    gene = entrez_ids,
    organism = "rno",
    pvalueCutoff = 0.05
  )
  
  # Guardar resultados
  write.csv(as.data.frame(ego), paste0("Homólogos/GO/", modelo, "_GO_results.csv"), row.names = FALSE)
  write.csv(as.data.frame(ekegg), paste0("Homólogos/KEGG/", modelo, "_KEGG_results.csv"), row.names = FALSE)
  
  # Guardar gráficos
  tiff(paste0("Homólogos/GO/", modelo, "_GO_dotplot.tiff"), width = 8, height = 6, units = "in", res = 600)
  print(dotplot(ego, showCategory = 20,font.size=8) + ggtitle(paste("GO -", modelo)))
  dev.off()
  
  tiff(paste0("Homólogos/KEGG/", modelo, "_KEGG_dotplot.tiff"), width = 8, height = 6, units = "in", res = 600)
  print(dotplot(ekegg, showCategory = 20,font.size=8) + ggtitle(paste("KEGG -", modelo)))
  dev.off()
  
  cat("Finalizado", modelo, "(", length(entrez_ids), "genes mapeados)\n\n")
}

# Ejecutar para ambos modelos
analisis_funcional("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/featureCounts/GSEA/DIO/DESeq2/DESeq2_results_DEGs.csv", "DIO")
analisis_funcional("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/featureCounts/GSEA/Zucker/DESeq2/DESeq2_results_DEGs.csv", "Zucker")
