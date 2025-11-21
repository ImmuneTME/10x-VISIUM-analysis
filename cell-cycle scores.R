library(SpaCET)
library(GiottoUtils)
library(Giotto)
library(GiottoClass)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(Seurat)
library(reticulate)

#load the giotto object with all (12 ChIO) samples preprocessed
results_folder = "S:/VISIUM/basic_expression/plots/"

py_config()
my_python_path = 'C:/Users/angel/Documents/.virtualenvs/r-reticulate/Scripts/python.exe'

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = TRUE,
                                  python_path = my_python_path)

#directories to the preprocessed giotto objects
#"S:/VISIUM/TFM/allsamples_giotto"
#"S:/VISIUM/TFM/allsamples15_giotto"

all_samples <- loadGiotto("S:/VISIUM/TFM/allsamples15_giotto",
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
#12 samples
s_allsamples_results[, c(7,8, 10:16, 54:56)] -> filtered_sresults
#15 samples
s_allsamples_results[, c(7,8, 15, 19:21)] -> filtered_sresults

#for the 12 samples with the 4 areas of interest:
# Rename columns 4 to 50 by removing the first 4 characters
colnames(filtered_sresults)[1:9] <- sub("^.{4}", "", colnames(filtered_sresults)[1:9])
filtered_sresults <- filtered_sresults %>%
  mutate(`area` = case_when(
    mTLS == 1 ~ "mTLSs",
    eTLS == 1 ~ "non-mTLSs",
    noTLS_B == 1 ~ "Dispersed B cells",
    `elses` == 1 ~ "Low B cell areas"
  ))

write.csv(filtered_sresults, "S:/VISIUM/cell-cycle/cycle_12samples.csv", row.names = TRUE, col.names = TRUE)
write.csv(filtered_sresults, "S:/VISIUM/cell-cycle/cycle_15samples.csv", row.names = TRUE, col.names = TRUE)
