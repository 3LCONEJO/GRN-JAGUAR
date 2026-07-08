# GRN Monocytes Project

This project organizes the GRN and Cicero/CellOracle workflow for monocytes.

## Project Map

1. `1_scripts/` - analysis scripts used to prepare inputs and run Cicero.
2. `2_data/` - raw and prepared data inputs, such as Seurat objects, GRN tables, and peaks.
3. `3_output/` - generated tables and analysis results, such as Cicero co-accessibility, CCANs, and Cytoscape exports.
4. `4_figures/` - plots and exported figures.

## Main Files

- `GRN.Rproj` - RStudio project file.
- `test_celloracle.ipynb` - Python notebook used for the CellOracle step.
- `1_scripts/get_seurat.r` - extracts peaks from the Seurat object.
- `1_scripts/create_cicero_cds.r` - builds the Cicero workflow outputs.
