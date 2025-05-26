# Cargar librerías
library(biomaRt)
library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Crear carpetas de salida
setwd("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/")
dir.create("Homologos", showWarnings = FALSE)
dir.create("Homologos/GO", showWarnings = FALSE)
dir.create("Homologos/KEGG", showWarnings = FALSE)

# Leer archivos de entrada
human_file <- "G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/Human_Datasets/Liver/GSEA/DESeq2/DESeq2_results_DEGs_NAFL.csv"
dio_file <- "G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/featureCounts/GSEA/DIO/DESeq2/DESeq2_results_DEGs.csv"
zucker_file <- "G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/featureCounts/GSEA/Zucker/DESeq2/DESeq2_results_DEGs.csv"

human_df <- read_csv(human_file)
dio_df <- read_csv(dio_file)
zucker_df <- read_csv(zucker_file)

# Filtrar genes significativos
human_sig <- subset(human_df, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)
dio_sig <- subset(dio_df, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)
zucker_sig <- subset(zucker_df, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)

# Extraer identificadores
human_ids <- unique(human_sig[,1])
dio_genes <- unique(dio_sig[,1])
zucker_genes <- unique(zucker_sig[,1])

# Conectar a Ensembl
#human_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
#rat_mart <- useEnsembl(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl")

#human_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
#rat_mart <- useEnsembl(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl", host = "useast.ensembl.org")


human_mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 105)
rat_mart <- useEnsembl(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl", version = 105)

# Obtener homólogos humano-rata
human_to_rat <- getLDS(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = human_ids,
  mart = human_mart,
  attributesL = c("external_gene_name", "ensembl_gene_id"),
  martL = rat_mart
)

colnames(human_to_rat) <- c("Human_Ensembl_ID", "Human_Symbol", "Rat_Gene_Symbol", "Rat_Ensembl_ID")

dio_genes<-unlist(dio_genes)
zucker_genes<-unlist(zucker_genes)

# Filtrar los que están en modelos DIO y Zucker
dio_homologs <- human_to_rat %>% filter(Rat_Gene_Symbol %in% dio_genes)
zucker_homologs <- human_to_rat %>% filter(Rat_Gene_Symbol %in% zucker_genes)

# Genes comunes entre los tres modelos
#genes_comunes <- intersect(dio_homologs$Rat_Gene_Symbol, zucker_homologs$Rat_Gene_Symbol)
#homologs_comunes <- human_to_rat %>% filter(Rat_Gene_Symbol %in% genes_comunes)

# Genes comunes entre modelos

homologs_comunes_DIO <- human_to_rat %>% filter(Rat_Gene_Symbol %in% dio_genes)
homologs_comunes_Zucker <- human_to_rat %>% filter(Rat_Gene_Symbol %in% zucker_genes)

# Genes humanos únicos para GO y KEGG
#genes_humanos_para_GO <- homologs_comunes %>%
# select(Human_Ensembl_ID, Human_Symbol) %>%
#distinct()



# Guardar resultados en CSV
write_csv(dio_homologs, "Homologos/dio_human_homologs_NAFL.csv")
write_csv(zucker_homologs, "Homologos/zucker_human_homologs_NAFL.csv")
#write_csv(homologs_comunes, "Homologos/genes_comunes_tres_modelos.csv")
#write_csv(genes_humanos_para_GO, "Homologos/genes_humanos_para_GO_KEGG.csv")



# Vías en DIO

# Genes humanos únicos para GO y KEGG
genes_humanos_para_GO_DIO <- homologs_comunes_DIO %>%
  dplyr::select(Human_Ensembl_ID, Human_Symbol) %>%
  distinct()

write_csv(genes_humanos_para_GO_DIO, "Homologos/genes_humanos_para_GO_KEGG_DIO_NAFL.csv")

# Convertir Ensembl a ENTREZ ID
entrez_ids <- bitr(
  genes_humanos_para_GO_DIO$Human_Ensembl_ID,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


# Análisis GO
go_results <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

write_csv(as.data.frame(go_results), "Homologos/GO/enrichGO_BP_results_DIO_NAFL.csv")

# Gráficos GO
tiff("Homologos/GO/dotplot_GO_DIO_NAFL.tiff", width = 8, height = 6, units = "in", res = 600)
dotplot(go_results, showCategory = 20,font.size=5)
dev.off()

tiff("Homologos/GO/barplot_GO_DIO_NAFL.tiff", width = 10, height = 6, units = "in", res = 600)
barplot(go_results, showCategory = 20)
dev.off()

if (nrow(go_results@result) > 0) {
  tiff("Homologos/GO/emapplot_GO_DIO_NAFL.tiff", width = 10, height = 8, units = "in", res = 600)
  emapplot(pairwise_termsim(go_results), showCategory = 20)
  dev.off()
  
  tiff("Homologos/GO/cnetplot_GO_DIO_NAFL.tiff", width = 12, height = 8, units = "in", res = 600)
  cnetplot(go_results, showCategory = 10, foldChange = NULL)
  dev.off()
}

# Análisis KEGG
kegg_results <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

kegg_readable <- setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_csv(as.data.frame(kegg_readable), "Homologos/KEGG/enrichKEGG_results_DIO_NAFL.csv")

# Gráficos KEGG
tiff("Homologos/KEGG/dotplot_KEGG_DIO_NAFL.tiff", width = 8, height = 6, units = "in", res = 600)
dotplot(kegg_readable, showCategory = 20)
dev.off()

tiff("Homologos/KEGG/barplot_KEGG_DIO_NAFL.tiff", width = 10, height = 6, units = "in", res = 600)
barplot(kegg_readable, showCategory = 20)
dev.off()


# Vías en Zucker

# Genes humanos únicos para GO y KEGG
genes_humanos_para_GO_Zucker <- homologs_comunes_Zucker %>%
  dplyr::select(Human_Ensembl_ID, Human_Symbol) %>%
  distinct()

write_csv(genes_humanos_para_GO_Zucker, "Homologos/genes_humanos_para_GO_KEGG_Zucker_NAFL.csv")

# Convertir Ensembl a ENTREZ ID
entrez_ids <- bitr(
  genes_humanos_para_GO_Zucker$Human_Ensembl_ID,
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)


# Análisis GO
go_results <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

write_csv(as.data.frame(go_results), "Homologos/GO/enrichGO_BP_results_Zucker_NAFL.csv")

# Gráficos GO
tiff("Homologos/GO/dotplot_GO_Zucker_NAFL.tiff", width = 8, height = 6, units = "in", res = 600)
dotplot(go_results, showCategory = 20, font.size=5)
dev.off()

tiff("Homologos/GO/barplot_GO_Zucker_NAFL.tiff", width = 10, height = 6, units = "in", res = 600)
barplot(go_results, showCategory = 20)
dev.off()

if (nrow(go_results@result) > 0) {
  tiff("Homologos/GO/emapplot_GO_Zucker_NAFL.tiff", width = 10, height = 8, units = "in", res = 600)
  emapplot(pairwise_termsim(go_results), showCategory = 20)
  dev.off()
  
  tiff("Homologos/GO/cnetplot_GO_Zucker_NAFL.tiff", width = 12, height = 8, units = "in", res = 600)
  cnetplot(go_results, showCategory = 10, foldChange = NULL)
  dev.off()
}

# Análisis KEGG
kegg_results <- enrichKEGG(
  gene = entrez_ids$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

kegg_readable <- setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
write_csv(as.data.frame(kegg_readable), "Homologos/KEGG/enrichKEGG_results_Zucker_NAFL.csv")

# Gráficos KEGG
tiff("Homologos/KEGG/dotplot_KEGG_Zucker_NAFL.tiff", width = 8, height = 6, units = "in", res = 600)
dotplot(kegg_readable, showCategory = 20)
dev.off()

tiff("Homologos/KEGG/barplot_KEGG_Zucker_NAFL.tiff", width = 10, height = 6, units = "in", res = 600)
barplot(kegg_readable, showCategory = 20)
dev.off()
