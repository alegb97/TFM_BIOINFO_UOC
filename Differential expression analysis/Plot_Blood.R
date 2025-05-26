library(readxl)
library(dplyr)
library(ggplot2)

# Crear carpetas de salida
setwd("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/Human_Datasets/Blood/")
output_dir <- "G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/Homologos/Blood/Plot"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cargar datos
counts <- as.data.frame(read_xlsx("GSE278204.xlsx", col_names = TRUE))
homologos_DIO <- read.csv("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/Homologos/dio_human_homologs_NAFL.csv")
homologos_Zucker <- read.csv("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/Homologos/zucker_human_homologs_NAFL.csv")
upregulated_DIO_Zucker <- read.csv("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/featureCounts/GSEA/OBS_Comparison_RRHO2/edgeR/upregulated_DIO_Zucker.csv")
downregulated_DIO_Zucker <- read.csv("G:/Mi unidad/Máster Bioinformática y Bioestadística/TFM/featureCounts/GSEA/OBS_Comparison_RRHO2/edgeR/downregulated_DIO_Zucker.csv")

# Pasar a mayúsculas
upregulated_DIO_Zucker <- toupper(upregulated_DIO_Zucker[,1])
downregulated_DIO_Zucker <- toupper(downregulated_DIO_Zucker[,1])

# Filtrar genes regulados por modelo correctamente
upregulated_genes_DIO <- homologos_DIO %>% filter(Human_Symbol %in% upregulated_DIO_Zucker)
downregulated_genes_DIO <- homologos_DIO %>% filter(Human_Symbol %in% downregulated_DIO_Zucker)

upregulated_genes_Zucker <- homologos_Zucker %>% filter(Human_Symbol %in% upregulated_DIO_Zucker)
downregulated_genes_Zucker <- homologos_Zucker %>% filter(Human_Symbol %in% downregulated_DIO_Zucker)

# Limpiar duplicados y calcular medias
counts$Gene <- as.character(counts$Gene)
counts <- counts %>%
  group_by(Gene) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = 'drop') %>%
  as.data.frame()
rownames(counts) <- counts[,1]

# Filtrar columnas de muestras
samples <- c(paste0("Control", 1:15), paste0("O-T2DM", 1:15))
counts <- counts[, c("Gene", samples)]

# Medias por grupo
mean_control <- rowMeans(counts[, 2:16])
mean_obese <- rowMeans(counts[, 17:31])

# Dataframe base para hacer gráfico
plot_df <- data.frame(
  Gene = counts$Gene,
  Mean_Control = mean_control,
  Mean_Obese = mean_obese
)

# Función  para hacer gráfico graficar
plot_model <- function(data, genes_ortologos, up_genes, down_genes, modelo_nombre) {
  data$Gene_upper <- toupper(data$Gene)
  genes_ortologos_upper <- toupper(genes_ortologos)
  up_genes_upper <- toupper(up_genes)
  down_genes_upper <- toupper(down_genes)
  
  data$Group <- "Other"
  data$Group[data$Gene_upper %in% genes_ortologos_upper] <- "Ortholog"
  data$Group[data$Gene_upper %in% down_genes_upper] <- "Downregulated Ortholog"
  data$Group[data$Gene_upper %in% up_genes_upper] <- "Upregulated Ortholog"
  
  data$Group <- factor(data$Group, levels = c("Other", "Ortholog", "Upregulated Ortholog", "Downregulated Ortholog"))
  
  data$Label <- ifelse(data$Group %in% c("Upregulated Ortholog", "Downregulated Ortholog"), data$Gene, NA)
  
  colores <- c(
    "Other" = "gray80",
    "Ortholog" = "orange",
    "Upregulated Ortholog" = "red",
    "Downregulated Ortholog" = "blue"
  )
  
  data$PointSize <- ifelse(data$Group %in% c("Upregulated Ortholog", "Downregulated Ortholog"), 4,
                           ifelse(data$Group == "Ortholog", 2.5, 1.5))
  
  # Evitar ceros para escala log
  data$Mean_Control <- pmax(data$Mean_Control, 0.01)
  data$Mean_Obese <- pmax(data$Mean_Obese, 0.01)
  
  p <- ggplot(data, aes(x = Mean_Control, y = Mean_Obese)) +
    geom_point(data = subset(data, Group == "Other"),
               aes(color = Group, size = PointSize),
               alpha = 0.6) +
    geom_point(data = subset(data, Group == "Ortholog"),
               aes(color = Group, size = PointSize),
               alpha = 0.85) +
    geom_point(data = subset(data, Group == "Upregulated Ortholog"),
               aes(color = Group, size = PointSize),
               alpha = 0.85) +
    geom_point(data = subset(data, Group == "Downregulated Ortholog"),
               aes(color = Group, size = PointSize),
               alpha = 0.85) +
    scale_size_identity() +
    geom_text(aes(label = Label),
              size = 4,
              hjust = 0.5,
              vjust = -0.7,
              na.rm = TRUE,
              fontface = "bold",
              check_overlap = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = colores, name = "Orthologs") +
    scale_x_log10(limits = c(0.01, 1000)) +
    scale_y_log10(limits = c(0.01, 1000)) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Human Blood Expression -", modelo_nombre),
      x = "Mean FPKM (Control)",
      y = "Mean FPKM (Obese)"
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
  
  
  ortologos_presentes <- data %>%
    filter(Group != "Other") %>%
    select(Gene, Mean_Control, Mean_Obese, Group)
  
  write.csv(ortologos_presentes, 
            file = file.path(output_dir, paste0("orthologs_in_plot_", modelo_nombre, ".csv")),
            row.names = FALSE)
  
  ggsave(filename = file.path(output_dir, paste0("scatterplot_", modelo_nombre, ".png")),
         plot = p, width = 9, height = 7, dpi = 300)
  
  return(p)
}

# Generar gráficos con genes regulados correctamente asignados
p1 <- plot_model(plot_df,homologos_DIO$Human_Symbol, upregulated_genes_DIO$Human_Symbol, downregulated_genes_DIO$Human_Symbol, "DIO")
p2 <- plot_model(plot_df,homologos_Zucker$Human_Symbol, upregulated_genes_Zucker$Human_Symbol, downregulated_genes_Zucker$Human_Symbol, "Zucker")

print(p1)
print(p2)
