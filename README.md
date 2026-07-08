# GRN Monocytes Project

This project organizes the GRN and Cicero/CellOracle workflow for monocytes.

## Project Map

1. `1_scripts/` - analysis scripts used to prepare inputs and run Cicero.
2. `2_data/` - raw and prepared data inputs, such as Seurat objects, GRN tables, and peaks.
3. `3_output/` - generated tables and analysis results, such as Cicero co-accessibility, CCANs, and Cytoscape exports.
4. `4_figures/` - plots and exported figures.

## Overview

The workflow follows this order:

1. Prepare peaks from the Seurat ATAC object.
2. Build the Cicero CDS and calculate co-accessibility links.
3. Export CCANs, gene activity matrices, and Cytoscape-ready networks.
4. Save plots and final figures in the figure folder.

## Main Files

- `GRN.Rproj` - RStudio project file.
- `test_celloracle.ipynb` - Python notebook used for the CellOracle step.
- `1_scripts/get_seurat.r` - extracts peaks from the Seurat object.
- `1_scripts/create_cicero_cds.r` - builds the Cicero workflow outputs.

## Notes

- Keep new analysis outputs in `3_output/`.
- Keep publication figures and exported plots in `4_figures/`.
- Use the numbered folders to keep the project easy to scan and share.