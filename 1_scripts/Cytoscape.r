# Cargar las librerías necesarias del tidyverse
suppressPackageStartupMessages({
  library(readr)
  library(tidyr)
  library(dplyr)
  library(RCy3)
})

print("Iniciando conversión a formato Cytoscape en R...")
csv_path <- "~/JAGUAR/GRN/data_folder/Base_GRN_Monocitos_JAGUAR.csv"

if (file.exists(csv_path)) {
  print("Cargando matriz CSV...")
  base_GRN <- read_csv(csv_path, show_col_types = FALSE)
  
  # 2. Transformación y FILTRADO TOP N
  print("Aplastando la matriz y extrayendo el TOP 5000...")
  cytoscape_edges <- base_GRN %>%
    pivot_longer(
      cols = -c(peak_id, gene_short_name, coaccess),
      names_to = "TF",
      values_to = "has_motif"
    ) %>%
    filter(has_motif == 1) %>%
    select(source = TF, target = gene_short_name, weight = coaccess) %>%
    # Eliminar flechas duplicadas del mismo TF al mismo Gen
    distinct(source, target, .keep_all = TRUE) %>%
    # ORDENAR de mayor a menor fuerza física (coaccess)
    arrange(desc(weight)) %>%
    # EXTRAER SOLO LAS MEJORES 5000 CONEXIONES
    slice_head(n = 5000) %>%
    mutate(interaction = "regulates")
  
  print(paste("Total de flechas únicas (controlado):", nrow(cytoscape_edges)))
  
  # --- NUEVO: Calcular atributos de los nodos (Conexiones y Tipo) ---
  print("Calculando jerarquía de nodos (Hubs)...")
  
  # Contamos cuántas conexiones salientes tiene cada TF
  origenes <- cytoscape_edges %>%
    group_by(source) %>%
    summarise(conexiones = n()) %>%
    rename(id = source) %>%
    mutate(tipo = "TF")
  
  # Los genes diana que no son TFs tendrán 1 conexión para el cálculo de tamaño
  destinos <- data.frame(id = setdiff(cytoscape_edges$target, origenes$id),
                         conexiones = 1,
                         tipo = "Diana")
  
  # Unimos ambos en la tabla final de nodos
  cytoscape_nodes <- bind_rows(origenes, destinos)
  # ------------------------------------------------------------------
  
  # 3. Enviar a Cytoscape
  print("========================================================")
  print("Conectando con Cytoscape...")
  
  tryCatch({
    cytoscapePing()
    
    createNetworkFromDataFrames(nodes = cytoscape_nodes, 
                                edges = cytoscape_edges, 
                                title = "Red_Monocitos_Top5000", 
                                collection = "Proyecto_JAGUAR")
    
    print("¡Red enviada! Creando Estilo Visual JAGUAR...")
    
    # --- NUEVO: Estilo Visual Programático ---
    style_name <- "Estilo_JAGUAR"
    createVisualStyle(style_name)
    setVisualStyle(style_name)
    
    # Forma circular sin bordes negros
    setNodeShapeDefault("ELLIPSE", style.name = style_name)
    setNodeBorderWidthDefault(0, style.name = style_name)
    
    # Color: Naranja (#FF5722) para TFs, Gris claro (#E0E0E0) para Dianas
    setNodeColorMapping('tipo', c('TF', 'Diana'), c('#FF5722', '#E0E0E0'), mapping.type = 'd', style.name = style_name)
    
    # Tamaño de Nodos basado en conexiones (Hubs gigantes, dianas pequeñas)
    min_con <- min(cytoscape_nodes$conexiones)
    max_con <- max(cytoscape_nodes$conexiones)
    setNodeSizeMapping('conexiones', c(min_con, max_con), c(15, 120), mapping.type = 'c', style.name = style_name)
    
    # Tamaño del Texto: Oculta el texto de genes menores (tamaño 1), resalta los TFs (tamaño 40)
    setNodeFontSizeMapping('conexiones', c(min_con, max_con), c(1, 40), mapping.type = 'c', style.name = style_name)
    setNodeLabelColorDefault('#000000', style.name = style_name)
    
    # Aristas: Gris muy claro y sin punta de flecha para no ensuciar el gráfico
    setEdgeColorDefault('#EEEEEE', style.name = style_name)
    setEdgeLineWidthDefault(0.5, style.name = style_name)
    setEdgeTargetArrowShapeDefault('NONE', style.name = style_name)
    # -----------------------------------------
    
    print("Aplicando diseño Prefuse Force-Directed (ideal para hubs)...")
    layoutNetwork('prefuse-force-directed')
    
    # 4. Guardar imagen
    image_path <- file.path(path.expand("~/JAGUAR/GRN/data_folder"), "Red_Monocitos_Top5000")
    exportImage(image_path, type = "PNG", resolution = 300)
    print(paste("¡Imagen guardada exitosamente en:", image_path, ".png !"))
    
  }, error = function(e) {
    print("ERROR con Cytoscape. Detalle:")
    print(e$message)
  })
  
} else {
  print("ERROR: Archivo no encontrado.")
}