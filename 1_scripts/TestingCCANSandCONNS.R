# Script: TestingCCANSandCONNS.R
# Objective: Try to merge the information of 
#            co-accessibility and cis-co-accessibility networks

# SECTION 0:
# Install Libraries
# ==============================================================================
# if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
# 
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# 
# if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
# 
# if (!requireNamespace("SeuratWrappers", quietly = TRUE)) pak::pak("satijalab/seurat-wrappers")
# 
# if (!requireNamespace("Signac", quietly = TRUE)) pak::pak("dgrunwald/signac") 
# 
# if (!requireNamespace("monocle3", quietly = TRUE)) pak::pak("cole-trapnell-lab/monocle3")
# 
# if (!requireNamespace("cicero", quietly = TRUE)) devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3", force = TRUE)
# 
# BiocManager::install(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats",
#                        "limma", "lme4", "S4Vectors", "SingleCellExperiment",
#                        "SummarizedExperiment", "batchelor", "HDF5Array",
#                        "ggrastr"), ask = FALSE)
# SECTION 1:
# Load Libraries and Define Paths
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(SeuratWrappers)
  library(monocle3)
  library(cicero) 
  library(EnsDb.Hsapiens.v86)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(VariantAnnotation)
  library(igraph)
  library(dplyr)
  library(ggforce)
  library(ggplot2)
})

project_root <- file.path(Sys.getenv("HOME"), "JAGUAR/GRN")
output_dir <- file.path(project_root, "3_output/cicero_mono_m3")
figure_dir <- file.path(project_root, "4_figures")
conns <- read.csv(file.path(output_dir, "monocitos_coaccess_cicero.csv"))
ccans <- read.csv(file.path(output_dir, "monocitos_CCANs.csv"))

# SECTION 2:
# Descriptive statistics
# ==============================================================================

print(head(ccans))

cat("||||How Many CCANs Exists?||||")
cat("There are",length(unique(ccans$CCAN)),"CCANs")

cat("||||How Many Peaks per CCANs Exists?||||")
hist(ccans$CCAN, breaks = length(unique(ccans$CCAN)))


empty <- c()

ccans_ids <- unique(ccans$CCAN)

region_df <- data.frame()

# This Step creates a data frame per CCAN where only kept the start and end of the ccan
# instead of having all the peaks.

for(i in ccans_ids){
  
  sub <- subset(ccans, CCAN == i)
  
  coords <- do.call(
    rbind,
    strsplit(sub$Peak, "-")
  )
  
  region_df <- rbind(
    region_df,
    data.frame(
      ccan = i,
      chromosome = unique(coords[,1]),
      start = min(as.numeric(coords[,2])),
      end = max(as.numeric(coords[,3])),
      n_peaks = nrow(sub)
    )
  )
}

region_df$size_bp <- region_df$end - region_df$start
region_df$size_kb <- region_df$size_bp / 1000
region_df$size_mb <- region_df$size_bp / 1e6

# Do bigger CCANs contain more peaks?

plot(
  region_df$size_kb,
  region_df$n_peaks,
  xlab = "Size (kb)",
  ylab = "N Peaks"
)

# Answer: Yes, Like it should?

# Top 20 ccans by number of peaks (picudos)
head(
  region_df[order(-region_df$n_peaks), ],
  20
)

# Top 20 ccans by size (grandotes)
head(
  region_df[order(-region_df$size_bp), ],
  20
)

# Section 3:
# Annotate ccans
# ==============================================================================

gr_ccan <- GRanges(
  seqnames = region_df$chromosome,
  ranges = IRanges(
    start = region_df$start,
    end = region_df$end
  )
)
gene_gr <- genes(EnsDb.Hsapiens.v86)

seqlevelsStyle(gene_gr) <- "UCSC"

hits <- findOverlaps(
  gr_ccan,
  gene_gr
)

ccan_genes <- data.frame(
  ccan = region_df$ccan[queryHits(hits)],
  gene = gene_gr$gene_name[subjectHits(hits)]
)

# Test to find OAS family genes

subset(
  ccan_genes,
  gene %in% c("OAS1","OAS2","OAS3","OASL")
)

# see how many genes are per ccan

genes_per_ccan <- aggregate(
  gene ~ ccan,
  ccan_genes,
  function(x) length(unique(x))
)

head(
  genes_per_ccan[
    order(-genes_per_ccan$gene),
  ],
  20
)

# Subset the ccan that contain oas

oas_peaks <- subset(
  ccans,
  CCAN == 1025
)$Peak

oas_conns <- subset(
  conns,
  Peak1 %in% oas_peaks &
    Peak2 %in% oas_peaks
)

# Create a graph of the connexions

g <- graph_from_data_frame(
  oas_conns[, c("Peak1","Peak2")],
  directed = FALSE
)

deg <- degree(g)

hubs <- names(sort(deg, decreasing = TRUE)[1:20])

# Make a data frame that contain the center of each peak

hub_df <- do.call(
  rbind,
  lapply(hubs, function(x){
    
    y <- strsplit(x, "-")[[1]]
    
    data.frame(
      chr = y[1],
      start = as.numeric(y[2]),
      end = as.numeric(y[3]),
      center = (as.numeric(y[2]) + as.numeric(y[3]))/2
    )
  })
)

# Given the overall coaccess per peak sum all that and obtain the total strength of the node

node_strength <- bind_rows(
  oas_conns %>%
    select(Peak = Peak1, coaccess),
  
  oas_conns %>%
    select(Peak = Peak2, coaccess)
) %>%
  group_by(Peak) %>%
  summarise(
    degree = n(),
    mean_coaccess = mean(coaccess),
    total_strength = sum(coaccess)
  )

node_strength |> arrange(desc(total_strength))

top10 <- node_strength |>
  arrange(desc(total_strength)) |>
  head(10)

peak_gr_peak_top <- GRanges(
  "chr12",
  IRanges(
    113040427,
    113041384
  )
)

findOverlaps(
  peak_gr_peak_top,
  gene_gr
)

dists <- distanceToNearest(
  peak_gr_peak_top,
  gene_gr,
  select = "all"
)

peaks_split <- do.call(rbind, strsplit(top10$Peak, "-"))

hubs_gr <- GRanges(
  seqnames = peaks_split[,1],
  ranges = IRanges(start = as.numeric(peaks_split[,2]), 
                   end = as.numeric(peaks_split[,3]))
)

top10_nearest_hits <- distanceToNearest(
  hubs_gr,
  gene_gr
)

oas_hits <- distanceToNearest(
  hubs_gr,
  subset(
    gene_gr,
    gene_name %in% c("OAS1","OAS2","OAS3")
  ),
  select = "all"
)

oas_conns |>
  arrange(desc(coaccess)) |>
  head(20)

nearest_gene <- gene_gr$gene_name[
  subjectHits(top10_nearest_hits)
]

nearest_distance <- mcols(top10_nearest_hits)$distance

summary_table <- data.frame(
  Peak = top10$Peak,
  Degree = top10$degree,
  Mean_Coaccess = round(top10$mean_coaccess,3),
  Total_Strength = round(top10$total_strength,3),
  Nearest_Gene = nearest_gene,
  Distance = nearest_distance
)

summary_table


coords <- do.call(
  rbind,
  strsplit(summary_table$Peak, "-")
)

summary_table$chr <- coords[,1]
summary_table$start <- as.numeric(coords[,2])
summary_table$end <- as.numeric(coords[,3])

oas_gr <- subset(
  gene_gr,
  gene_name %in% c("OAS1","OAS2","OAS3")
)

hits_oas <- overlapsAny(
  hubs_gr,
  oas_gr
)

summary_table$Overlaps_OAS <- hits_oas

oas_locus <- GRanges(
  "chr12",
  IRanges(
    start = min(start(oas_gr)),
    end = max(end(oas_gr))
  )
)

dist_oas <- distanceToNearest(
  hubs_gr,
  oas_locus
)

summary_table$Distance_to_OAS <- mcols(dist_oas)$distance

summary_table$HubScore <-
  summary_table$Degree *
  summary_table$Mean_Coaccess


# Section 4
# Filter coaccess by significance
# ==============================================================================

threshold <- quantile(
  oas_conns$coaccess,
  0.95
)

oas_sig <- subset(
  oas_conns,
  coaccess >= threshold
)

oas_sig <- oas_sig[oas_sig$Peak1 < oas_sig$Peak2, ]

peak_center <- function(peak) {
  
  x <- strsplit(peak, "-")[[1]]
  
  start <- as.numeric(x[2])
  end   <- as.numeric(x[3])
  
  (start + end) / 2
}

arc_df <- oas_sig

arc_df$x1 <- sapply(
  arc_df$Peak1,
  peak_center
)

arc_df$x2 <- sapply(
  arc_df$Peak2,
  peak_center
)

arc_df$curvature <- rescale(
  arc_df$coaccess,
  to = c(0.1, 0.8)
)

arc_df$alpha_val <- rescale(arc_df$coaccess, to = c(0.45, 0.95))
arc_df$width_val <- rescale(arc_df$coaccess, to = c(0.8, 2.8))

# SECTION 5
# Plot coaccessibility + peaks + genes
# ==============================================================================
# Panel 1 (top)  -> Arch
# Panel 2 (mid)   -> Peaks
# Panel 3 (below)   -> Genes

# --------------------------------------------------------------------------
# 0. Whole Window Size 
# --------------------------------------------------------------------------

ccan_region <- subset(region_df, ccan == ccan_id)

region_chr   <- ccan_region$chromosome
region_start <- ccan_region$start - 5000
region_end   <- ccan_region$end + 5000

# --------------------------------------------------------------------------
# 1. Peak Panel Plot
# --------------------------------------------------------------------------

peaks_split_all <- do.call(rbind, strsplit(oas_peaks, "-"))

peaks_df <- data.frame(
  chr   = peaks_split_all[, 1],
  start = as.numeric(peaks_split_all[, 2]),
  end   = as.numeric(peaks_split_all[, 3])
)

peaks_df$width  <- peaks_df$end - peaks_df$start
peaks_df$is_hub <- oas_peaks %in% top10$Peak

peak_plot <- ggplot(peaks_df) +
  geom_rect(
    aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.8, fill = is_hub)
  ) +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), guide = "none") +
  coord_cartesian(xlim = c(region_start, region_end)) +
  labs(title = "ATAC peaks (red = top10 hubs)") +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)))

# --------------------------------------------------------------------------
# 2. Arch Plot Panel
# --------------------------------------------------------------------------

arc_plot <- ggplot() +
  coord_cartesian(xlim = c(region_start, region_end), ylim = c(0, 1))

for (i in seq_len(nrow(arc_df))) {
  
  arc_plot <- arc_plot +
    geom_curve(
      data = arc_df[i, ],
      aes(x = x1, y = 0, xend = x2, yend = 0),
      curvature  = -arc_df$curvature[i],
      alpha      = arc_df$alpha_val[i],
      linewidth  = arc_df$width_val[i],
      color      = "steelblue4"
    )
}

arc_plot <- arc_plot +
  labs(title = paste0("Co-Accsessibility (>= cuantil 0.95, coaccess >= ",
                      round(threshold, 3), ")")) +
  theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)))

# --------------------------------------------------------------------------
# 3. Gene Plot Panel
# --------------------------------------------------------------------------

genes_in_region <- subset(
  gene_gr,
  as.character(seqnames(gene_gr)) == region_chr &
    start(gene_gr) < region_end &
    end(gene_gr)   > region_start
)

genes_df <- as.data.frame(genes_in_region)
genes_df$is_oas <- genes_df$gene_name %in% c("OAS1", "OAS2", "OAS3")

genes_df$gene_name <- factor(
  genes_df$gene_name,
  levels = unique(
    genes_df$gene_name[
      order(genes_df$start)
    ]
  )
)

genes_plus  <- subset(genes_df, strand == "+")
genes_minus <- subset(genes_df, strand == "-")

gene_plot <- ggplot() +
  { if (nrow(genes_plus) > 0)
    geom_segment(
      data = genes_plus,
      aes(x = start, xend = end, y = gene_name, yend = gene_name, color = is_oas),
      linewidth = 3,
      arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "open")
    )
  } +
  { if (nrow(genes_minus) > 0)
    geom_segment(
      data = genes_minus,
      aes(x = start, xend = end, y = gene_name, yend = gene_name, color = is_oas),
      linewidth = 3,
      arrow = arrow(length = unit(0.2, "cm"), ends = "first", type = "closed")
    )
  } +
  scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "firebrick"), guide = "none") +
  coord_cartesian(xlim = c(region_start, region_end)) +
  labs(
    title = "Genes in Region",
    x = paste0("Position: ", region_chr, " (bp)")
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)),
    axis.title.y = element_blank()
  )

# --------------------------------------------------------------------------
# 4. Combine Plots
# --------------------------------------------------------------------------

final_plot <- arc_plot / peak_plot / gene_plot +
  plot_layout(heights = c(2, 0.6, 1.8))

final_plot

ggsave(
  filename = file.path(figure_dir, paste0("CCAN_", ccan_id, "_track_plot.png")),
  plot = final_plot,
  width = 9, height = 6.5, dpi = 300
)

# SECTION 3 & 4 & 5:
# Generalized Annotation, Filtering, and Plotting
# ==============================================================================
library(patchwork) # For combining plots
library(scales)    # For rescale()

# 1. Prepare Genomic Ranges for Genes
gene_gr <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gene_gr) <- "UCSC"

# 2. Define the Generalized Plotting Function
plot_ccan_track <- function(ccan_id, ccans_df, conns_df, region_df, gene_gr, fig_dir) {
  
  # --- Subsetting Data for Current CCAN ---
  current_peaks <- subset(ccans_df, CCAN == ccan_id)$Peak
  current_conns <- subset(conns_df, Peak1 %in% current_peaks & Peak2 %in% current_peaks)
  
  # Safety check: If no connections exist, skip plotting
  if(nrow(current_conns) < 1) {
    message(paste("Skipping CCAN", ccan_id, "- No connections found."))
    return(NULL)
  }
  
  # --- Calculate Hubs ---
  # Total strength of the node
  node_strength <- bind_rows(
    current_conns %>% select(Peak = Peak1, coaccess),
    current_conns %>% select(Peak = Peak2, coaccess)
  ) %>%
    group_by(Peak) %>%
    summarise(
      degree = n(),
      mean_coaccess = mean(coaccess),
      total_strength = sum(coaccess)
    )
  
  top10 <- node_strength %>% arrange(desc(total_strength)) %>% head(10)
  
  # Extract hubs GRanges to find nearest genes (to highlight them dynamically)
  peaks_split <- do.call(rbind, strsplit(top10$Peak, "-"))
  hubs_gr <- GRanges(
    seqnames = peaks_split[,1],
    ranges = IRanges(start = as.numeric(peaks_split[,2]), 
                     end = as.numeric(peaks_split[,3]))
  )
  
  # Find genes nearest to the top 10 hubs to highlight them in the plot
  nearest_hits <- distanceToNearest(hubs_gr, gene_gr)
  hub_genes_to_highlight <- unique(gene_gr$gene_name[subjectHits(nearest_hits)])
  
  # --- Filter Connections by Significance ---
  threshold <- quantile(current_conns$coaccess, 0.95, na.rm = TRUE)
  current_sig <- subset(current_conns, coaccess >= threshold)
  
  # Safety check: If threshold is too high and removes everything, fallback to top connections
  if(nrow(current_sig) == 0) current_sig <- head(arrange(current_conns, desc(coaccess)), 10)
  
  current_sig <- current_sig[current_sig$Peak1 < current_sig$Peak2, ]
  
  peak_center <- function(peak) {
    x <- strsplit(peak, "-")[[1]]
    (as.numeric(x[2]) + as.numeric(x[3])) / 2
  }
  
  arc_df <- current_sig
  arc_df$x1 <- sapply(arc_df$Peak1, peak_center)
  arc_df$x2 <- sapply(arc_df$Peak2, peak_center)
  
  # Handle scaling safely if there's only 1 row
  if(nrow(arc_df) > 1) {
    arc_df$curvature <- rescale(arc_df$coaccess, to = c(0.1, 0.8))
    arc_df$alpha_val <- rescale(arc_df$coaccess, to = c(0.45, 0.95))
    arc_df$width_val <- rescale(arc_df$coaccess, to = c(0.8, 2.8))
  } else {
    arc_df$curvature <- 0.45; arc_df$alpha_val <- 0.7; arc_df$width_val <- 1.5
  }
  
  # --- Prepare Genomic Window ---
  ccan_region <- subset(region_df, ccan == ccan_id)
  if(nrow(ccan_region) == 0) return(NULL)
  
  region_chr   <- ccan_region$chromosome
  chr_folder <- file.path(fig_dir, region_chr)
  if(!dir.exists(chr_folder)) dir.create(chr_folder, recursive = TRUE)
  
  region_start <- ccan_region$start - 5000
  region_end   <- ccan_region$end + 5000
  
  # --- PLOT 1: Peaks ---
  peaks_split_all <- do.call(rbind, strsplit(current_peaks, "-"))
  peaks_df <- data.frame(
    chr   = peaks_split_all[, 1],
    start = as.numeric(peaks_split_all[, 2]),
    end   = as.numeric(peaks_split_all[, 3])
  )
  peaks_df$is_hub <- current_peaks %in% top10$Peak
  
  peak_plot <- ggplot(peaks_df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.8, fill = is_hub)) +
    scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), guide = "none") +
    coord_cartesian(xlim = c(region_start, region_end)) +
    labs(title = paste("ATAC peaks (red = top hubs for CCAN", ccan_id, ")")) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)))
  
  # --- PLOT 2: Arcs ---
  arc_plot <- ggplot() + coord_cartesian(xlim = c(region_start, region_end), ylim = c(0, 1))
  
  for (i in seq_len(nrow(arc_df))) {
    arc_plot <- arc_plot +
      geom_curve(
        data = arc_df[i, ],
        aes(x = x1, y = 0, xend = x2, yend = 0),
        curvature  = -arc_df$curvature[i],
        alpha      = arc_df$alpha_val[i],
        linewidth  = arc_df$width_val[i],
        color      = "steelblue4"
      )
  }
  
  arc_plot <- arc_plot +
    labs(title = paste0("CCAN ", ccan_id, " Co-Accessibility (>= 0.95 quant: ", round(threshold, 3), ")")) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)))
  
  # --- PLOT 3: Genes ---
  genes_in_region <- subset(
    gene_gr,
    as.character(seqnames(gene_gr)) == region_chr &
      start(gene_gr) < region_end &
      end(gene_gr)   > region_start
  )
  
  genes_df <- as.data.frame(genes_in_region)
  
  if(nrow(genes_df) > 0) {
    genes_df$is_hub_gene <- genes_df$gene_name %in% hub_genes_to_highlight
    genes_df$gene_name <- factor(genes_df$gene_name, levels = unique(genes_df$gene_name[order(genes_df$start)]))
    genes_plus  <- subset(genes_df, strand == "+")
    genes_minus <- subset(genes_df, strand == "-")
    
    gene_plot <- ggplot() +
      { if (nrow(genes_plus) > 0)
        geom_segment(data = genes_plus, aes(x = start, xend = end, y = gene_name, yend = gene_name, color = is_hub_gene),
                     linewidth = 3, arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "open"))
      } +
      { if (nrow(genes_minus) > 0)
        geom_segment(data = genes_minus, aes(x = start, xend = end, y = gene_name, yend = gene_name, color = is_hub_gene),
                     linewidth = 3, arrow = arrow(length = unit(0.2, "cm"), ends = "first", type = "closed"))
      } +
      scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "firebrick"), guide = "none") +
      coord_cartesian(xlim = c(region_start, region_end)) +
      labs(title = "Genes in Region (red = nearest to hubs)", x = paste0("Position: ", region_chr, " (bp)")) +
      theme_minimal(base_size = 9) +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 9, hjust = 0, margin = margin(b = 2)),
            axis.title.y = element_blank())
  } else {
    # Empty placeholder if no genes fall in this window
    gene_plot <- ggplot() + theme_void() + labs(title = "No genes mapped in this region")
  }
  
  # --- COMBINE AND SAVE ---
  
  num_genes <- nrow(genes_df)
  
  dynamic_height <- 5 + (num_genes * 0.25)
  dynamic_height <- max(7, dynamic_height)
  
  gene_panel_weight <- max(1.8, num_genes * 0.12)
  
  final_plot <- arc_plot / peak_plot / gene_plot + 
    plot_layout(heights = c(2, 0.6, gene_panel_weight))
  
  file_name <- file.path(chr_folder, paste0("CCAN_", ccan_id, "_track_plot.pdf"))
  # Asegúrate de definir bien la ruta donde están tus 2328 PDFs ahora mismo
figure_dir <- file.path(project_root, "4_figures") 

for (i in 1:nrow(region_df)) {
  
  # Extraer la información de la fila actual
  ccan_id <- region_df$ccan[i]
  chrom <- region_df$chromosome[i]
  
  # Definir cómo se llama el archivo y dónde debería ir
  file_name <- paste0("CCAN_", ccan_id, "_track_plot.pdf")
  origen <- file.path(figure_dir, file_name)
  destino_dir <- file.path(figure_dir, chrom)
  destino_final <- file.path(destino_dir, file_name)
  
  # Lógica de movimiento
  if (file.exists(origen)) {
    
    # Crear la carpeta del cromosoma si aún no existe
    if (!dir.exists(destino_dir)) {
      dir.create(destino_dir, recursive = TRUE)
    }
    
    # Mover el archivo
    file.rename(from = origen, to = destino_final)
  }
}
  ggsave(
    filename = file_name, 
    plot = final_plot, 
    width = 11, 
    height = dynamic_height, 
    limitsize = FALSE
  )
  
  message(paste("Saved PDF plot for CCAN:", ccan_id, "| Genes:", num_genes, "| Height:", round(dynamic_height, 1)))
}


# ==============================================================================
# SECTION 6: Execute Loop Across All CCANs
# ==============================================================================

# Create figure directory if it doesn't exist
if(!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)

all_ccan_ids <- unique(ccans$CCAN)
total_ccans <- length(all_ccan_ids)

cat("Starting to generate plots for", total_ccans, "CCANs...\n")

for (i in 1:total_ccans) {
  current_id <- all_ccan_ids[i]
  
  # Wrap in tryCatch so one bad CCAN doesn't break the whole loop
  tryCatch({
    plot_ccan_track(
      ccan_id   = current_id, 
      ccans_df  = ccans, 
      conns_df  = conns, 
      region_df = region_df, 
      gene_gr   = gene_gr, 
      fig_dir   = figure_dir
    )
  }, error = function(e) {
    message(paste("Error processing CCAN", current_id, ":", e$message))
  })
}

cat("Finished processing all CCANs.\n")