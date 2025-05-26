# Cargar librerías necesarias
library(Rsubread)

# Obtener lista de archivos BAM ordenados
bam_files <- list.files(pattern = "DIO_S[0-9]+_sorted.bam$")

# Ejecutar featureCounts para DIO
conteo_DIO <- featureCounts(files = bam_files,
                            annot.ext = "rn6.refGene.gtf",
                            isGTFAnnotationFile = TRUE,
                            isPairedEnd = TRUE,
                            strandSpecific = 2,  #  Reversely stranded 
                            nthreads = 12)

# Guardar resultados en CSV
write.csv(conteo_DIO$counts, "featureCounts_DIO_results.csv")
write.csv(conteo_DIO$stat, "featureCounts_DIO_stats.csv")