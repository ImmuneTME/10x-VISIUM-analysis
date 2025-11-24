library(SpaCET)
library(GiottoUtils)
library(Giotto)
library(GiottoClass)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(Seurat)
library(reticulate)

#load the giotto object with all  samples preprocessed
results_folder = "S://"

py_config()
my_python_path = 'C://python.exe'

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = TRUE,
                                  python_path = my_python_path)


all_samples <- loadGiotto("C:/",
                          load_params = list(),
                          reconnect_giottoImage = TRUE,
                          python_path = my_python_path,
                          init_gobject = TRUE,
                          verbose = TRUE)

#transform that giotto object into a seurat object
giottoToSeuratV5(
  gobject = all_samples,
  tech = "Visium",
  res_type = "lowres",
  verbose = TRUE,
  set_defaults = TRUE,
  simplify = TRUE
) -> s_all_samples

#a list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat
#extract the cell-cycle gene markers and compute scoring
s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes

s_all_samples <- CellCycleScoring(s_all_samples, 
                               s.features = s.genes, 
                               g2m.features = g2m.genes, 
                               set.ident = TRUE) 


s_all_samples@meta.data -> s_allsamples_results
#filter the output table. column indexes differ by giotto input
s_allsamples_results[, c()] -> filtered_sresults


write.csv(filtered_sresults, "C:/.csv", row.names = TRUE, col.names = TRUE)


