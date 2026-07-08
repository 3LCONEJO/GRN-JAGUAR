
mono_seurat <- readRDS("~/JAGUAR/GRN/2_data/mono_seurat_filtered.rds")

# 1. Extraer los nombres de las filas (que son los picos de consenso)
picos_crudos <- rownames(mono_seurat[["ATAC"]])

# 2. Convertir el formato de Signac (chr1-12345-67890) al formato de CellOracle (chr1_12345_67890)
# Reemplazamos los guiones por guiones bajos
picos_formateados <- gsub("-", "_", picos_crudos)

# 3. Crear el DataFrame estructurado
picos_df <- data.frame(peak_id = picos_formateados)

# 4. Guardar como CSV para enviarlo directamente a tu Jupyter Notebook
write.csv(picos_df, "~/JAGUAR/GRN/2_data/peaks.csv", row.names = FALSE, quote = FALSE)

print(paste("[DEBUG] -> Se exportaron exitosamente", nrow(picos_df), "picos unificados."))