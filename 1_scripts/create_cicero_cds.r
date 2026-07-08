#==============================================================================
# SCRIPT PURPOSE:
# This script creates a Cicero Cell Data Set (CDS) object from a Seurat
# object containing single-cell ATAC-seq data (specifically for monocytes).
# Cicero is used to predict cis-regulatory interactions (e.g.,
# enhancer-promoter links) by analyzing co-accessibility patterns.
#==============================================================================

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) pak::pak("satijalab/seurat-wrappers")

if (!requireNamespace("Signac", quietly = TRUE)) pak::pak("dgrunwald/signac") 

if (!requireNamespace("monocle3", quietly = TRUE)) pak::pak("cole-trapnell-lab/monocle3")

if (!requireNamespace("cicero", quietly = TRUE)) devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3", force = TRUE)

BiocManager::install(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats",
                       "limma", "lme4", "S4Vectors", "SingleCellExperiment",
                       "SummarizedExperiment", "batchelor", "HDF5Array",
                       "ggrastr"), ask = FALSE)

if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

#==============================================================================
# SECTION 1: LOAD LIBRARIES
#------------------------------------------------------------------------------
# Load all necessary packages for the analysis.
#==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratWrappers)
  library(monocle3)
  library(cicero) 
  library(EnsDb.Hsapiens.v86)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
})

#==============================================================================
# SECTION 1.5: DEFINE FILE PATHS
#------------------------------------------------------------------------------
project_root <- file.path(Sys.getenv("HOME"), "JAGUAR/GRN")
output_dir <- file.path(project_root, "3_output/cicero_mono_m3")
figure_dir <- file.path(project_root, "4_figures")
seurat_object_path <- file.path(project_root, "2_data/mono_seurat_filtered.rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
mono_seurat <- readRDS(seurat_object_path)

#==============================================================================
# SECTION 2: CREATE OR LOAD CICERO CDS OBJECT
#------------------------------------------------------------------------------
# This section checks if a Cicero CDS object has already been created and
# saved. If it exists, it loads it to save time. Otherwise, it generates it
# from the source Seurat object.
#==============================================================================

cicero_cds_path <- file.path(project_root, "3_output/cicero_mono")

# Check if the output directory with the saved Cicero object already exists.
if (dir.exists(cicero_cds_path)) {
  # If it exists, load the previously saved Cicero CDS object.
  print("[DEBUG] -> [Sección 2] Objeto Cicero encontrado, cargando desde disco...")
  cicero_cds <- load_monocle_objects(directory_path = cicero_cds_path)
  
} else {
  # If it does not exist, proceed with creating it from scratch.
  print("[DEBUG] -> [Sección 2] Objeto Cicero NO encontrado, creando uno nuevo...")
  
  # Step 2.1: Load Seurat
  mono_seurat <- readRDS("~/JAGUAR/GRN/2_data/mono_seurat_filtered.rds")
  
  # Step 2.2: Convert the Seurat object to a Monocle3 cell_data_set (CDS).
  mono_cds <- as.cell_data_set(x = mono_seurat)
  mono_cds <- cluster_cells(mono_cds, reduction_method = "UMAP")
  
  # Step 2.3: Extract the UMAP coordinates from the CDS object.
  umap_coords <- reducedDims(mono_cds)$UMAP
  
  # Step 2.4: Create the Cicero CDS object.
  cicero_cds <- make_cicero_cds(mono_cds,
                                reduced_coordinates = umap_coords, k = 50)
  
  # Step 2.5: Save the newly created Cicero CDS object to disk.
  print("[DEBUG] -> [Sección 2] Guardando objeto Cicero CDS en disco...")
  save_monocle_objects(cicero_cds, directory_path = cicero_cds_path)
}

# 1. Load the human genome reference (hg38) to provide genomic
# coordinates for Cicero.
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
genome_df <- data.frame("chr" = names(genome), "length" = genome)

# 2. Calculate co-accessibility scores using Cicero.
print("[DEBUG] -> [Sección 3] Calculando conexiones de co-accesibilidad...")
conns <- run_cicero(cicero_cds, genome_df, sample_num = 100) # Agregado sample_num=100 para estabilidad

# 3. Export the result for CellOracle
conns_filtered <- subset(conns, coaccess > 0)
write.csv(conns_filtered, file.path(output_dir, "monocitos_coaccess_cicero.csv"), row.names = FALSE, quote = FALSE)
print("[DEBUG] -> [Sección 3] Conexiones exportadas exitosamente.")

#==============================================================================
# SECTION 4: VISUALIZATION (HUMAN hg38 ADAPTATION VIA ENSDB)
#==============================================================================
print("[DEBUG] -> [Sección 4] Obteniendo anotaciones de genes humanos desde EnsDb.Hsapiens.v86...")

conns <- read.csv(file.path(output_dir, "monocitos_coaccess_cicero.csv"))

# Extract Annotation Info from EnsDb
tx <- transcripts(EnsDb.Hsapiens.v86, columns = c("gene_id", "tx_id", "gene_name"))
gene_anno <- as.data.frame(tx)

# Match column names
gene_anno$chromosome <- paste0("chr", gene_anno$seqnames)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$tx_id
gene_anno$symbol <- gene_anno$gene_name
print(paste("[DEBUG]   -> Anotaciones procesadas:", nrow(gene_anno), "transcritos."))

print("[DEBUG] -> [Sección 4] Generando gráfico de conexiones en PDF...")
pdf(file.path(figure_dir, "Cicero_Connections_OAS1_Locus.pdf"), width = 40, height = 20)
plot_connections(conns, "chr12", 112490000 , 113200000,
                 gene_model = gene_anno,
                 coaccess_cutoff = 0.1,
                 connection_width = 0.5,
                 alpha_by_coaccess = TRUE,
                 collapseTranscripts = "longest")
dev.off()
print("[DEBUG]   -> Gráfico guardado: Cicero_Connections_OAS1_Locus.pdf")

#==============================================================================
# SECTION 5: FINDING CIS-CO-ACCESSIBILITY NETWORKS (CCANs)
#==============================================================================
print("[DEBUG] -> [Sección 5] Buscando módulos CCANs (Comunidades Louvain)...")
CCAN_assigns <- generate_ccans(conns_filtered)
print(paste("[DEBUG]   -> Se asignaron CCANs a", nrow(CCAN_assigns), "picos."))
write.csv(CCAN_assigns, file.path(output_dir, "monocitos_CCANs.csv"), row.names = FALSE)

#==============================================================================
# SECTION 6: CICERO GENE ACTIVITY SCORES
#==============================================================================
print("[DEBUG] -> [Sección 6] Calculando 'Gene Activity Scores' (Traducción ATAC -> RNA)...")

# [Step 6.1] -> Prepare annotation with promoters
print("[DEBUG]   -> Extrayendo coordenadas de sitios de inicio de transcripción (TSS)...")
pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
pos <- pos[!duplicated(pos$transcript),] 
pos$end <- pos$start + 1 

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
neg <- neg[!duplicated(neg$transcript),] 
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)
gene_annotation_sub <- gene_annotation_sub[, c("chromosome", "start", "end", "symbol")]
names(gene_annotation_sub)[4] <- "gene"

# [Step 6.2] -> Annotate CDS with promoters
print("[DEBUG]   -> Anotando el objeto CDS con los TSS locales...")
mono_cds <- annotate_cds_by_site(mono_cds, gene_annotation_sub)

# [Step 6.3] -> Build matrix (not-normalized)
print("[DEBUG]   -> Construyendo matriz de actividad génica (Unnormalized)...")
unnorm_ga <- build_gene_activity_matrix(mono_cds, conns)
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, !Matrix::colSums(unnorm_ga) == 0]

# [Step 6.4] -> Bring metadata
print("[DEBUG]   -> Asegurando métrica de 'sitios totales accesibles por célula'...")
if (is.null(pData(mono_cds)$num_genes_expressed)) {
  if (!is.null(pData(mono_cds)$nFeature_ATAC)) {
    num_genes <- pData(mono_cds)$nFeature_ATAC
  } else {
    num_genes <- Matrix::colSums(counts(mono_cds) > 0)
  }
} else {
  num_genes <- pData(mono_cds)$num_genes_expressed
}
names(num_genes) <- row.names(pData(mono_cds))

# [Step 6.5] -> Normalization
print("[DEBUG]   -> Normalizando Gene Activity Scores...")
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)

print(paste("[DEBUG] -> [ÉXITO] Dimensiones finales de Actividad Génica (Genes x Células):", 
            paste(dim(cicero_gene_activities), collapse = " x ")))

# Save sparce matrix
saveRDS(cicero_gene_activities, file.path(output_dir, "cicero_gene_activities_matrix.rds"))

print("[DEBUG] ==============================================================")
print("[DEBUG] ¡TODOS LOS StepS DEL TUTORIAL FUERON COMPLETADOS CON ÉXITO!")
print("[DEBUG] ==============================================================")

#==============================================================================
# SECTION 7: INTEGRATE GENE ACTIVITY SCORES BACK INTO SEURAT
#==============================================================================
print("[DEBUG] -> [Sección 7] Integrando 'Gene Activity Scores' en el objeto Seurat...")

# 1. Load the saved gene activity matrix.
cicero_gene_activities <- readRDS(file.path(output_dir, "cicero_gene_activities_matrix.rds"))

# 2. Filter the matrix to ensure it exactly matches the cells in the Seurat object.
cicero_gene_activities <- cicero_gene_activities[, colnames(mono_seurat)]

# 3. Create a new Assay in the Seurat object to store the activity scores.
mono_seurat[["ACTIVITY"]] <- CreateAssayObject(counts = cicero_gene_activities)

# 4. Set the new "ACTIVITY" assay as the default for visualization.
DefaultAssay(mono_seurat) <- "ACTIVITY"

# 5. Visualize gene activity on the UMAP.
print("[DEBUG]   -> Generando FeaturePlot para genes de ejemplo...")
FeaturePlot(mono_seurat, features = c("MNT", "TTC7A", "CARS2"), pt.size = 0.5)

# 6. Example of searching for specific genes in the activity matrix.
coincidencias <- grep("CDC", rownames(cicero_gene_activities), ignore.case = TRUE, value = TRUE) 
print(coincidencias)

print("[DEBUG] -> [Sección 7] Integración y visualización completadas.")
