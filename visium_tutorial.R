library(Giotto)
library(umap)
library(scran)
library(MAST)
library(multinet)
library(RTriangle)
library(reticulate)
library(harmony)
library(tidyverse)
library(R.utils)
library(tiff)
library(biomaRt)
library(FactoMineR)
library(MCPcounter)
library(openxlsx)
library(quadprog)
library(Rfast)
library(Matrix)
library(spacexr)
#library(STdeconvolve)
#library(mistyR)
#install_miniconda()



#options(timeout = 6000000000) ### set this to avoid timeout error
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

#devtools::install_github("edsgard/trendsceek")
setwd("S:/VISIUM/basic_expression/")


library(reticulate)

use_condaenv("r-reticulate", required = TRUE)

reticulate::conda_list()
#use_python('C:/Users/angel/AppData/Local/r-miniconda/envs/r-reticulate/python.exe', required = TRUE)
#use_miniconda(condaenv = 'C:/Users/angel/AppData/Local/r-miniconda/envs/r-reticulate/python.exe', required = TRUE)

use_miniconda(condaenv = "r-reticulate", required = TRUE)



#options(timeout = 6000000000) ### set this to avoid timeout error
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

#devtools::install_github("edsgard/trendsceek")
setwd("S:/VISIUM/basic_expression/")

#1.- DATA PREPARATION
results_folder = "S:/VISIUM/basic_expression/plots/"

#py_config()
#my_python_path = 'C:/Users/angel/AppData/Local/r-miniconda/envs/r-reticulate/python.exe'
my_python_path = 'C:/Users/angel/Documents/.virtualenvs/r-reticulate/Scripts/python.exe'
#my_python_path = "C:/Users/lungc/miniconda3/python.exe" #python directory workstation
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = TRUE,
                                  python_path = my_python_path)

######## WARNING, IMAGES v2: Y AXIS FLIPPED. ELIMINATE THE "2" IF FUNCTIONALITY RETURNS

#1.1.- open each of the samples
#ChIO SAMPLES
results_dir_40 <- "S:/VISIUM/visium_results/outs_p40/"
visium_40 = createGiottoVisiumObject(visium_dir = results_dir_40,
                                     expr_data = 'filter',
                                     png_name = 'tissue_lowres_image.png',
                                     gene_column_index = 2,
                                     instructions = instrs)

results_dir_27 <- "S:/VISIUM/visium_results/outs_p27/"
visium_27 = createGiottoVisiumObject(visium_dir = results_dir_27,
                                     expr_data = 'filter',
                                     png_name = 'tissue_lowres_image.png',
                                     gene_column_index = 2,
                                     instructions = instrs)

results_dir_28 <- "S:/VISIUM/visium_results/outs_p28/"
visium_28 = createGiottoVisiumObject(visium_dir = results_dir_28,
                                     expr_data = 'filter',
                                     png_name = 'tissue_lowres_image.png',
                                     gene_column_index = 2,
                                     instructions = instrs)

results_dir_23 <- "S:/VISIUM/visium_results/outs_p23/"
visium_23 = createGiottoVisiumObject(visium_dir = results_dir_23,
                                     expr_data = 'filter',
                                     png_name = 'tissue_lowres_image.png',
                                     gene_column_index = 2,
                                     instructions = instrs)

results_dir_14 <- "S:/VISIUM/visium_results/Visium_N1_14_output/outs_p14/"
visium_14  = createGiottoVisiumObject(visium_dir = results_dir_14,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

results_dir_15 <- "S:/VISIUM/visium_results/Visium_N1_15_output/outs_p15/"
visium_15  = createGiottoVisiumObject(visium_dir = results_dir_15,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

results_dir_21 <- "S:/VISIUM/visium_results/Visium_N1_21_output/outs_p21/"
visium_21  = createGiottoVisiumObject(visium_dir = results_dir_21,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

results_dir_24 <- "S:/VISIUM/visium_results/Visium_N1_24_output/outs_p24/"
visium_24  = createGiottoVisiumObject(visium_dir = results_dir_24,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

results_dir_33 <- "S:/VISIUM/visium_results/Visium_N1_33_output/outs_p33/"
visium_33  = createGiottoVisiumObject(visium_dir = results_dir_33,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

results_dir_43 <- "S:/VISIUM/visium_results/Visium_N1_43_output/outs_p43/"
visium_43  = createGiottoVisiumObject(visium_dir = results_dir_43,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

results_dir_49 <- "S:/VISIUM/visium_results/Visium_N1_49_output/outs_p49/"
visium_49  = createGiottoVisiumObject(visium_dir = results_dir_49,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

results_dir_50 <- "S:/VISIUM/visium_results/Visium_N1_50_output/outs_p50/"
visium_50  = createGiottoVisiumObject(visium_dir = results_dir_50,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

results_dir_04 <- "S:/VISIUM/visium_results/Visium_N2_04_output/outs_p04/"
visium_04  = createGiottoVisiumObject(visium_dir = results_dir_04,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

results_dir_46 <- "S:/VISIUM/visium_results/Visium_N2_46_output/outs_p46/"
visium_46  = createGiottoVisiumObject(visium_dir = results_dir_46,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

results_dir_72 <- "S:/VISIUM/visium_results/Visium_N2_72_output/outs_p72/"
visium_72  = createGiottoVisiumObject(visium_dir = results_dir_72,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)

#1.2.- open and load the updated metadata onto each visium object
read.csv("S:/VISIUM/updated metadata/metadata_14.csv", header = TRUE) -> metadata_14
addCellMetadata(gobject = visium_14, new_metadata = metadata_14, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_14

read.csv("S:/VISIUM/updated metadata/metadata_15.csv", header = TRUE) -> metadata_15
addCellMetadata(gobject = visium_15, new_metadata = metadata_15, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_15

read.csv("S:/VISIUM/updated metadata/metadata_24.csv", header = TRUE) -> metadata_24
addCellMetadata(gobject = visium_24, new_metadata = metadata_24, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_24

read.csv("S:/VISIUM/updated metadata/metadata_40.csv", header = TRUE) -> metadata_40
addCellMetadata(gobject = visium_40, new_metadata = metadata_40, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_40

read.csv("S:/VISIUM/updated metadata/metadata_23.csv", header = TRUE) -> metadata_23
addCellMetadata(gobject = visium_23, new_metadata = metadata_23, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_23

read.csv("S:/VISIUM/updated metadata/metadata_43.csv", header = TRUE) -> metadata_43
addCellMetadata(gobject = visium_43, new_metadata = metadata_43, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_43

read.csv("S:/VISIUM/updated metadata/metadata_49.csv", header = TRUE) -> metadata_49
addCellMetadata(gobject = visium_49, new_metadata = metadata_49, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_49

read.csv("S:/VISIUM/updated metadata/metadata_50.csv", header = TRUE) -> metadata_50
addCellMetadata(gobject = visium_50, new_metadata = metadata_50, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_50

read.csv("S:/VISIUM/updated metadata/metadata_21.csv", header = TRUE) -> metadata_21
addCellMetadata(gobject = visium_21, new_metadata = metadata_21, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_21

read.csv("S:/VISIUM/updated metadata/metadata_27.csv", header = TRUE) -> metadata_27
addCellMetadata(gobject = visium_27, new_metadata = metadata_27, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_27

read.csv("S:/VISIUM/updated metadata/metadata_28.csv", header = TRUE) -> metadata_28
addCellMetadata(gobject = visium_28, new_metadata = metadata_28, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_28

read.csv("S:/VISIUM/updated metadata/metadata_33.csv", header = TRUE) -> metadata_33
addCellMetadata(gobject = visium_33, new_metadata = metadata_33, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_33

read.csv("S:/VISIUM/updated metadata/metadata_04.csv", header = TRUE) -> metadata_04
addCellMetadata(gobject = visium_04, new_metadata = metadata_04, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_04

read.csv("S:/VISIUM/updated metadata/metadata_46.csv", header = TRUE) -> metadata_46
addCellMetadata(gobject = visium_46, new_metadata = metadata_46, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_46

read.csv("S:/VISIUM/updated metadata/metadata_72.csv", header = TRUE) -> metadata_72
addCellMetadata(gobject = visium_72, new_metadata = metadata_72, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_72

all_samples <- joinGiottoObjects(
  gobject_list = list(visium_14, visium_15, visium_24, visium_40, visium_23, visium_43,
                      visium_49, visium_50, visium_21, visium_27, visium_28, visium_33, visium_04, visium_46, visium_72
  ),
  gobject_names = c("14", "15", "24", "40", "23", "43", "49", "50", "21", "27", "28", "33", "04", "46", "72"
  ),
  join_method = "shift"  #only one to keep spatial data separate among samples (padding)
)
 gc()

#2- QUALITY CONTROL
#2.1.- remove spots with zero genes detected 
all_samples <- filterGiotto(gobject = all_samples,
                            expression_threshold = 1,
                            feat_det_in_min_cells = 1, #min num of spots (cell) that need to express a gene (feature)
                            min_det_feats_per_cell = 1, #min num of genes (features) that need to be detected in a spot (cell)
                            expression_values = c('raw'),
                            verbose = T)
#2.2.- normalize the merged object
all_samples <- normalizeGiotto(gobject = all_samples, scalefactor = 10000, verbose = T)

#all_samples <- normalizeGiotto(gobject = all_samples, scalefactor = 10000, norm_methods = "standard", logbase = 2)

#2.3.- actual filter
all_samples <- filterGiotto(gobject = all_samples,
                            expression_threshold = 1,
                            feat_det_in_min_cells = 5, #min num of spots (cell) that need to express a gene (feature)
                            min_det_feats_per_cell = 300, #min num of genes (features) that need to be detected in a spot (cell)
                            expression_values = c('normalized'),
                            verbose = T)
gc()
#all_samples <- addStatistics(gobject = all_samples)

#norm_matrix <- all_samples@expression[["cell"]][["rna"]][["normalized"]]@exprMat
#norm_df     <- as.data.frame(as.matrix(norm_matrix))

#3.- Add deconvolution metadata
#read.csv("S:/VISIUM/ALL_MCP.csv", header = TRUE, row.names = 1) -> metadata_mcp
#addCellMetadata(gobject = all_samples, new_metadata = metadata_mcp, by_column = TRUE, column_cell_ID = "cell_ID") -> all_samples

pDataDT(all_samples) -> metadata_all

metadata_all <- metadata_all %>%
  mutate(elses = case_when(
    TLS == 0 & noTLS_B == 0 ~ 1,
    TLS == 0 & noTLS_B == 1 ~ 0,
    TLS == 1 ~ 0
  ))

table(metadata_all$elses)
addCellMetadata(gobject = all_samples, new_metadata = metadata_all, by_column = TRUE, column_cell_ID = "cell_ID") -> all_samples

read.csv("S:/VISIUM/deconvolution/analysis/RESULTS/no_malignant/SpaCET_nomalignant.csv", header = TRUE, row.names = 1) -> spacetdeconv
addCellMetadata(gobject = all_samples, new_metadata = spacetdeconv, by_column = TRUE, column_cell_ID = "cell_ID") -> all_samples
#write.csv(metadata_all, "processed_allsamples_metadata.csv", col.names = TRUE) 

rm(list = c("metadata_14", "metadata_15", "metadata_21", "metadata_23", "metadata_24", "metadata_27",
            "metadata_28", "metadata_33", "metadata_40", "metadata_43", "metadata_49", "metadata_50"))
rm(list = c("visium_14", "visium_15", "visium_21", "visium_23", "visium_24", "visium_27",
            "visium_28", "visium_33", "visium_40", "visium_43", "visium_49", "visium_50"))
gc()

#trendSceek(all_samples)

#4.- SUBSET
# Extract logical indices for the subset, then add its corresponding metadata and split it from the merged samples
TLS_indices <- all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$TLS ==  1 #| all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B == 1
TLS_metadata <- all_samples@cell_metadata$cell$rna@metaDT[TLS_indices, ]
TLS_metadata$cell_ID -> TLS_spots
TLS_samples <- subsetGiotto(all_samples, cell_ids = TLS_spots)

eTLS_indices <- all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$eTLS ==  1 #| all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B == 1
eTLS_metadata <- all_samples@cell_metadata$cell$rna@metaDT[eTLS_indices, ]
eTLS_metadata$cell_ID -> eTLS_spots
eTLS_samples <- subsetGiotto(all_samples, cell_ids = eTLS_spots)

mTLS_indices <- all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$mTLS ==  1 #| all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B == 1
mTLS_metadata <- all_samples@cell_metadata$cell$rna@metaDT[mTLS_indices, ]
mTLS_metadata$cell_ID -> mTLS_spots
mTLS_samples <- subsetGiotto(all_samples, cell_ids = mTLS_spots)

noTLS_indices <- all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$TLS ==  0 #| all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B == 1
noTLS_metadata <- all_samples@cell_metadata$cell$rna@metaDT[noTLS_indices, ]
noTLS_metadata$cell_ID -> noTLS_spots
noTLS_samples <- subsetGiotto(all_samples, cell_ids = noTLS_spots)

CPR_indices <- all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$response == "CPR" #| all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B == 1
CPR_metadata <- all_samples@cell_metadata$cell$rna@metaDT[CPR_indices, ]
CPR_metadata$cell_ID -> CPR_spots
CPR_samples <- subsetGiotto(all_samples, cell_ids = CPR_spots)

NCPR_indices <- all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$response == "NCPR" #| all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B == 1
NCPR_metadata <- all_samples@cell_metadata$cell$rna@metaDT[NCPR_indices, ]
NCPR_metadata$cell_ID -> NCPR_spots
NCPR_samples <- subsetGiotto(all_samples, cell_ids = NCPR_spots)

noTLS_B_indices <- all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B ==  1 #| all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B == 1
noTLS_B_metadata <- all_samples@cell_metadata$cell$rna@metaDT[noTLS_B_indices, ]
noTLS_B_metadata$cell_ID -> noTLS_B_spots
noTLS_B_samples <- subsetGiotto(all_samples, cell_ids = noTLS_B_spots)

else_indices <- else_indices <- all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$elses == 1 #| all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B == 1
else_metadata <- all_samples@cell_metadata$cell$rna@metaDT[else_indices, ]
else_metadata$cell_ID -> else_spots
else_samples <- subsetGiotto(all_samples, cell_ids = else_spots)

saveGiotto(
  gobject = mTLS_samples ,
  foldername = "mTLSs_giotto",
  dir = "S:/VISIUM/TFM/",
  method = c("RDS"),
  method_params = list(),
  overwrite = FALSE,
  export_image = TRUE,
  image_filetype = "PNG",
  include_feat_coord = TRUE,
  verbose = TRUE,
  metadata = TRUE
)

loadGiotto(
  path_to_folder = "S:/VISIUM/TFM/mTLSs_giotto",
  reconnect_giottoImage = TRUE,
  python_path = 'C:/Users/angel/Documents/.virtualenvs/r-reticulate/Scripts/python.exe',
  init_gobject = TRUE,
  verbose = TRUE
) -> mTLStest

#pDataDT(all_samples) -> regulators
#5- DIFFERENTIAL EXPRESSION ANALYSIS
#5.1 SCRAN

table(numbercheck$response)
table(eTLS_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["response"]])
table(else_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["list_ID"]])
diff_expr_results8 <- findScranMarkers(gobject = mTLS_samples,
                                      expression_values = "normalized",
                                      cluster_column = "response", 
                                      group_1 = "CPR", 
                                      group_1_name = "CPR",
                                      group_2 = "NCPR", 
                                      group_2_name = "NCPR")
diff_expr_results1[[1]] -> mTLS_POV
diff_expr_results1[[2]] -> other_POV
write.csv(noTLSB_POV, "S:/VISIUM/basic_expression/SCRAN/group vs all others/SCRAN_noTLSBvsothers_all_noTLSBpov.csv", col.names = TRUE)
write.csv(rest_POV, "S:/VISIUM/basic_expression/SCRAN/group vs all others/SCRAN_noTLSBvsothers_all_otherspov.csv", col.names = TRUE)
getwd()

#5.2. PLOTS
#GENES SELECTION

topgenes_scran = CPR_POV[, head(.SD, 30), by = 'cluster']$feats #takes top 70 DEGs

TLS <- c("MS4A1", "CD19", "CR2", "CXCL13")
HLA <- c("HHLA1", "HHLA2", "HHLA3", "HLA-A", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA",
         "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1",
         "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-E", "HLA-F", "HLA-G")
IG <- c("IGBP1", "IGDCC3", "IGDCC4", "IGF1", "IGF1R", "IGF2BP1", "IGF2BP2", "IGF2R", 
        "IGFALS", "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", 
        "IGFBPL1", "IGFL1", "IGFL2", "IGFL3", "IGFL4", "IGFLR1", "IGFN1", "IGHA1", 
        "IGHA2", "IGHD", "IGHE", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHJ2", "IGHJ6", 
        "IGHM", "IGHMBP2", "IGHV1-46", "IGHV1-58", "IGHV5-10-1", "IGHV6-1", "IGIP", 
        "IGKC", "IGKV1D-42", "IGKV4-1", "IGKV5-2", "IGKV6D-41", "IGLC1", "IGLC7",
        "IGLJ1", "IGLJ5", "IGLJ6", "IGLL1", "IGLV2-33", "IGLV3-1", "IGLV3-22", "IGLV3-32",
        "IGLV4-60", "IGLV6-57", "IGLV9-49", "IGSF1", "IGSF10", "IGSF11", "IGSF21", 
        "IGSF22", "IGSF23", "IGSF3", "IGSF5", "IGSF6", "IGSF8", "IGSF9", "IGSF9B")
cxcr5 <- c("CXCR5")

TLS_markers <- list(TLS_Van = c("CD3", "CD4", "CD8", "CD20", "MS4A1", "CD21", "CD23"), 
                    TLS_Cruz = c("MS4A1", "CD19", "CR2", "CXCL13"),
                    TLS_meylan = c("APOE", "C1QA", "C7", "CD52", "CD79A", "CXCL12", "DERL3", "FCRL5", "IGHA1", "IGHG1", 
                                   "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM", "IGKC", "IGLC1", "IGLC2", "IGLC3", 
                                   "IL7R", "JCHAIN", "LUM", "MZB1", "PIM2", "PTGDS", "PTLP", "SSR4", "TRBC2", "XBP1"),
                    TLS_mature = c("CD4", "PD1", "CXCR5", "CXCR13", "BCL6", "CD20", "MS4A1", "CD23", "PNAd"),
                    citokines = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")
)

Bsmarkers <- c("CXCR4", "MS4A1", "CD22", "CXCL13", "BCL2")
ATOCHA <- c("SMARCA4", "RB1", "DICER1")

involution <- c("MKI67", "BCL6")


#heatmaps
#plotMetaDataHeatmap(all_samples, selected_feats = topgenes_scran, metadata_cols = c('leiden_clus'))
#plotMetaDataHeatmap(TLS_samples, selected_feats = topgenes_scran, metadata_cols = c('combo'))
plotMetaDataHeatmap(TLS_samples, selected_feats = involution, metadata_cols = c('response'), expression_values = "normalized",
                    x_text_size = 16, y_text_size = 16, strip_text_size = 20 )

ggsave(
  filename = "heatmap_allsamples_TLS.tiff",
  plot = last_plot(),
  path = "S:/VISIUM/basic_expression/plots/",
  units = "in",
  width = 898 / 72,
  height = 1436 / 72, 
  dpi = 300 #increase to 300 if resolution needs to be corrected
)


violinPlot(TLS_samples, feats = involution, cluster_column = "response", strip_position = "top", 
           color_violin = "cluster", cluster_color_code = c("#0B645A", 'grey'), expression_values = "normalized",
           strip_text = 24, axis_text_x_size = 17, axis_text_y_size = 17)

'#0B645A'
"#A91E45FF"


ggsave(
  filename = "violin_all_response.tiff",
  plot = last_plot(),
  path = "S:/VISIUM/ATOCHA/",
  units = "in",
  width = 898 / 72,
  height = 2436 / 72, 
  dpi = 300 #increase to 300 if resolution needs to be corrected
)

ggsave(
  filename = "violin_NCPR_TLS.tiff",
  plot = last_plot(),
  path = "S:/VISIUM/ATOCHA/",
  units = "in",
  width = 898 / 72,
  height = 1436 / 72, 
  dpi = 300 #increase to 300 if resolution needs to be corrected
)

##### PLOTTING

read.csv("S:/VISIUM/basic_expression/SCRAN/TLS_SCRAN_CPRpov.csv", header = 1, row.names = 1) -> CPR_POV


#volcano plots
library(ggrepel)
volcano_data <- data.frame(gene = CPR_POV$feats, logFC = CPR_POV$summary.logFC, adjpval = -log10(CPR_POV$FDR), FDR = CPR_POV$FDR)
HLA_volcano_data <- subset(volcano_data, gene %in% HLA)
IG_volcano_data <- subset(volcano_data, gene %in% IG)

logFC_threshold <- 0.5
fdr_threshold <- 0.001

#TEXT fdr and logFC thresholds
#ggplot(volcano_data, aes(x = logFC, y = adjpval)) + 
#  geom_point(aes(color = FDR < fdr_threshold), alpha = 0.6) +
#  labs(x = "Log Fold Change", y = "-log10(p-value)") +
#  scale_color_manual(
#    name = paste("FDR < ", fdr_threshold),
#    values = c("grey", "red")) +
#  theme_minimal() +
#  geom_text_repel(data = subset(volcano_data, abs(logFC) > logFC_threshold & FDR <fdr_threshold),
#                  aes(label = gene), size = 4)

#TEXT logFC threshold only
ggplot(volcano_data, aes(x = logFC, y = adjpval)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.6) +
  labs(x = "Log Fold Change", y = "-log10(FDR)") +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() + 
  theme(
    axis.title.x = element_text(size = 16),  # Increase x-axis title size
    axis.title.y = element_text(size = 16),  # Increase y-axis title size
    legend.text = element_text(size = 12),   # Increase legend text size
    legend.title = element_text(size = 14)   # Increase legend title size
  )+ #TILL NOW WAS THE ACTUAL VOLCANO
  geom_text_repel(data = subset(volcano_data, abs(logFC) > logFC_threshold & FDR < 0.001),
                                    aes(label = gene), size = 4) #TAGS
ggsave(
  filename = "CPRPOV_noTLSB_HLA.tiff",
  plot = last_plot(),
  path = "C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/plots/FDR axis volcano/noTLS_B/",
  units = "in",
  width = 975 / 72,
  height = 1036 / 72, 
  dpi = 72 #increase to 300 if resolution needs to be corrected
)


rm(list = c("volcano_data", "NO_POV", "YES_POV", "diff_expr_results"))

##gene enrichment analysis
#open up DEG results (scran)
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/TLS_SCRAN_CPRpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/TLS_SCRAN_NOpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/eTLS_SCRAN_NOpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/eTLS_SCRAN_CPRpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/mTLS_SCRAN_CPRpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/mTLS_SCRAN_NOpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/bcells_SCRAN_NOpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/bcells_SCRAN_CPRpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/noTLS_B_SCRAN_CPRpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/noTLS_B_SCRAN_NOpov.csv"

#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/noTLSBvseTLS_SCRAN_eTLSpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/TLS_SCRAN_mTLSpov.csv"
#"C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/SCRAN/noTLSBvsTLS_SCRAN_TLSpov.csv"


read.csv("S:/VISIUM/basic_expression/SCRAN/else_SCRAN_CPRpov.csv", header = 1, row.names = 1) -> de_genes2

read.csv("S:/VISIUM/basic_expression/SCRAN/mTLS_SCRAN_NCPRpov.csv", header = 1, row.names = 1) -> de_genes

# Filter genes based on FDR or logFC 
#FILTER ORA ONLY
significant_genes <- de_genes[de_genes$FDR < 0.01 & de_genes$summary.logFC > 0, "feats"]
#all_genes <- de_genes[de_genes$FDR < 0.05, "feats"] #gsea
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
library(vctrs)
library(GOplot)
library(enrichplot)
library(DOSE)

# Define gene list (convert if necessary, e.g., to Entrez IDs)
ego <- enrichGO(gene = significant_genes, #top 1000 DEGs
                OrgDb = org.Hs.eg.db,   # For human data; replace for other species
                keyType = "SYMBOL", 
                ont = "BP",             # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                pAdjustMethod = "fdr", 
                pvalueCutoff = 0.05, 
                qvalueCutoff = 0.05)

options(enrichplot.colours = c('#DD3497', "#F3E96B"))

enrichplot::dotplot(ego, showCategory = 15) + 
  ggtitle("GO: BP overexpressed routes in non-progressor B cells") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),           # Increase title size
    axis.title.x = element_text(size = 14),                        # Increase x-axis title size
    axis.text.x = element_text(size = 12),                         # Increase x-axis value size
    axis.text.y = element_text(size = 16.5),                         # Adjust font size for route names
    legend.title = element_text(size = 12),                        # Increase legend title size
    legend.text = element_text(size = 10),                         # Increase legend text size
    plot.margin = unit(c(1, 1, 1, 3), "cm")                       # Increase left margin for longer names
  )

ggsave(
  filename = "GO_MF_NO_bcells.tiff",
  plot = last_plot(),
  path = "C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/plots/GO & GSEA/",
  units = "in",
  width = 1138 / 72,
  height = 1236 / 72, 
  dpi = 300 #increase to 300 if resolution needs to be corrected
)

barplot(ego, showCategory = 15) +
  ggtitle("GO: MF overexpressed routes in non-progressor B cells") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),           # Increase title size
    axis.title.x = element_text(size = 14),                        # Increase x-axis title size
    axis.text.x = element_text(size = 12),                         # Increase x-axis value size
    axis.text.y = element_text(size = 16.5),                         # Adjust font size for route names
    legend.title = element_text(size = 12),                        # Increase legend title size
    legend.text = element_text(size = 10),                         # Increase legend text size
    plot.margin = unit(c(1, 1, 1, 3), "cm")                       # Increase left margin for longer names
  )

ggsave(
  filename = "GO_MF_NO_bcells_bar.tiff",
  plot = last_plot(),
  path = "C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/plots/GO & GSEA/",
  units = "in",
  width = 1138 / 72,
  height = 1236 / 72, 
  dpi = 300 #increase to 300 if resolution needs to be corrected
)

#GSEA - USE ALL GENES AVAILABLE
gene_ranking <- de_genes$summary.logFC
names(gene_ranking) <- de_genes$feats
gene_ranking <- sort(gene_ranking, decreasing = TRUE)
head(gene_ranking) #quickcheck

gsea_results_go <- gseGO(geneList = gene_ranking,
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = "BP",           # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                         pAdjustMethod = "fdr",
                         pvalueCutoff = 0.05,
                         verbose = FALSE)

write.csv(gsea_results_go, "S:/VISIUM/basic_expression/GSEA/results/GSEA_BP_else_CPRpov.csv" )

#gsea_results_go@result -> gsea_df

#gsea plotting once its been calculated - revisiting results
read.csv("S:/VISIUM/basic_expression/SCRAN/TLS_SCRAN_CPRpov.csv", header = 1, row.names = 1) -> de_genes

gene_ranking <- de_genes$summary.logFC
names(gene_ranking) <- de_genes$feats
gene_ranking <- sort(gene_ranking, decreasing = TRUE)
head(gene_ranking) #quickcheck


read.csv("S:/VISIUM/basic_expression/GSEA/results/GSEA_BP_TLS_CPR.csv", header = TRUE) -> gsea

gsea_results_go <- new("gseaResult",
                       result = gsea,
                       geneList = gene_ranking,    # If you want to include the gene list, pass it here
                       organism = "Homo sapiens",
                       setType = "BP",
                       keytype = "SYMBOL",
                       readable = FALSE)

options(enrichplot.colours = c('#DD3497', "#F3E96B"))
#options(enrichplot.colours = c( "orange", "#1A75BA"))
#c('#FFF7F3', '#FDE0DD', '#FCC5C0', '#FA9FB5', '#F768A1', '#DD3497', '#AE017E', '#7A0177', '#49006A'))
selected_pathways_BP <- c("cell killing", "positive regulation of programmed cell death", "positive regulation of apoptotic process",
                          "antigen processing and presentation", "antigen processing and presentation of peptide antigen via MHC class II",
                          "lymphocyte proliferation", "lympchocyte differentiation", "positive regulation of lymphocyte activation",
                          "B cell mediated immunity", "B cell receptor signaling pathway",
                          "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",
                          "T cell proliferation", "T cell differentiation", "positive regulation of T cell activation",
                          "positive regulation of cytokine production")

dotplot(gsea_results_go, showCategory = selected_pathways_BP, x = "NES") +
  scale_size_area(max_size = 30) + #dot size keeping proportions/scale
  ggtitle("Upregulated Biological Processes in CPR mTLS") +
  theme(
    plot.title = element_text(size = 19, face = "bold"),           # Increase title size
    axis.title.x = element_text(size = 17),                        # Increase x-axis title size
    axis.text.x = element_text(size = 15),                         # Increase x-axis values size
    axis.text.y = element_text(size = 19),                         # Adjust font size for route names
    legend.title = element_text(size = 14),                        # Increase legend title size
    legend.text = element_text(size = 14),                         # Increase legend text size
    plot.margin = unit(c(1, 1, 1, 3), "cm")                       # Increase left margin for longer names
  ) 

ggsave(
  filename = "NES_selected_exclusive_BP_CPR_mTLS_dotplot.tiff",
  plot = last_plot(),
  path = "C:/Users/angel/Documents/VISIUM ANALYSIS/basic_expression/GSEA/plots/selected pathways/",
  units = "in",
  width = 1138 / 72,
  height = 1236 / 72,   
  dpi = 300 #increase to 300 if resolution needs to be corrected
)


#for(feat in TLS){ spatFeatPlot2D(gobject = all_samples, expression_values = "normalized", feats = feat, point_size = 2) }
####TISSUE / SAMPLE PLOTTING
##############################
alphabet = c('0' = '#72A2C0', '1' = '#F2917D', '2' = '#BC589B', '3' = '#7F7F7F', '4' = '#ABDDDE', '5' = '#FDD262', '6'='#EB545C','7'='#00475F', '8' = '#FFFF0A','9' = '#D7EF9B','10' = '#9E0142','11' = '#F2A104','12' = '#289E92','13' = '#80FF08')
colors1 = c("#192E5B","#1D65A6","#72A2C0","#F3E96B","#F2A104")
colors2 = c('#F7FCF0', '#E0F3DB', '#CCEBC5', '#A8DDB5', '#7BCCC4', '#4EB3D3', '#2B8CBE', '#0868AC', '#084081')
colors3 = c('#FFF7F3', '#FDE0DD', '#FCC5C0', '#FA9FB5', '#F768A1', '#DD3497', '#AE017E', '#7A0177', '#49006A')
custom_colors1 <- c("#B25D91FF", "#CB87B4FF", "#EFC7E6FF", "#1BB6AFFF", "#088BBEFF", "#172869FF")

custom_colors2 <- c("#EFC7E6FF", "#088BBEFF")
custom_colors3 <- c("lightyellow", "blue")


spots14 <- grep("^(14-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)
spots15 <- grep("^(15-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)
spots24 <- grep("^(24-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)
spots40 <- grep("^(40-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)

spots23 <- grep("^(23-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)
spots43 <- grep("^(43-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)
spots49 <- grep("^(49-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)
spots50 <- grep("^(50-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)

spots21 <- grep("^(21-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)
spots27 <- grep("^(27-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)
spots28 <- grep("^(28-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)
spots33 <- grep("^(33-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)

sub14 <- subsetGiottoLocs(all_samples, x_max = 23314, x_min = 4741, return_gobject = TRUE)
sub15 <- subsetGiottoLocs(all_samples, x_max = 50466, x_min = 29772, return_gobject = TRUE)
sub24 <- subsetGiottoLocs(all_samples, x_max = 75723, x_min = 56643, return_gobject = TRUE)
sub40 <- subsetGiottoLocs(all_samples, x_max = 102761, x_min = 83427, return_gobject = TRUE)

sub23 <- subsetGiottoLocs(all_samples, x_max = 130887, x_min = 110472, return_gobject = TRUE)
sub43 <- subsetGiottoLocs(all_samples, x_max = 158266, x_min = 138905, return_gobject = TRUE)
sub49 <- subsetGiottoLocs(all_samples, x_max = 184459, x_min = 164837, return_gobject = TRUE)
sub50 <- subsetGiottoLocs(all_samples, x_max = 211954, x_min = 192103, return_gobject = TRUE)

sub21 <- subsetGiottoLocs(all_samples, x_max = 239178, x_min = 219802, return_gobject = TRUE)
sub27 <- subsetGiottoLocs(all_samples, x_max = 266093, x_min = 245394, return_gobject = TRUE)
sub28 <- subsetGiottoLocs(all_samples, x_max = 293124, x_min = 272408, return_gobject = TRUE)
sub33 <- subsetGiottoLocs(all_samples, x_max = 319049, x_min = 299445, return_gobject = TRUE)

(names(all_samples@images)) -> image_names

for(feat in ATOCHA){ spatFeatPlot2D(gobject = sub24, show_image = T, image_name = "24-image",
                                expression_values = "normalized", feats = feat, point_size = 2) }

spatFeatPlot2D(gobject = sub33, show_image = T, image_name = "33-image", 
               expression_values = "normalized", feats = ATOCHA, point_size = 2,
               cow_n_col = 3, point_alpha = 0.8, cell_color_gradient = custom_colors3)


#METADATA
#all samples, all images
spatPlot2D(gobject = all_samples,  show_image = T, image_name = image_names, point_size = 2,
           point_alpha = 10, cell_color = "T.Cells", color_as_factor = F, 
           cell_color_gradient = custom_colors1)

spatPlot2D(gobject = all_samples,  show_image = T, image_name = "40-image", point_size = 4,
           point_alpha = 10, cell_color = "total_expr", color_as_factor = F, 
           select_cells = spots40, show_other_cells = FALSE,
           cell_color_gradient = custom_colors1)

#certain sample
spatPlot2D(gobject = sub33,  show_image = T, image_name = "33-image", point_size = 4,
           point_alpha = 10, cell_color = "total_expr", color_as_factor = F, 
           select_cells = spots33, show_other_cells = FALSE, coord_fix_ratio = 1,
           cell_color_gradient = custom_colors1)

spatPlot2D(gobject = sub40,  show_image = F, image_name = "40-image", point_size = 4,
           point_alpha = 10, cell_color = "B.lineage", color_as_factor = F, 
           select_cells = spots40, show_other_cells = FALSE, coord_fix_ratio = 1,
           cell_color_gradient = custom_colors3, background_color = "black")
#CELL-TYPE GRADIENT
spatPlot2D(gobject = sub43,  show_image = T, point_size = 4, image_name = "43-image",
           point_alpha = 10, cell_color = "GC", cell_color_gradient = custom_colors3,
           color_as_factor = FALSE, background_color = "black")

#CELL-TYPE BINARY
spatPlot2D(gobject = sub40,  show_image = T, image_name = "40-image", point_size = 4,
           point_alpha = 10, cell_color = "TLS", cell_color_gradient = colors1,
           color_as_factor = FALSE, background_color = "white")

#CELL-TYPE BINARY
spatPlot2D(gobject = sub40,  show_image = FALSE, point_size = 4,
           cell_color = "TLS", #cell_color_gradient = colors1, point_alpha = 10,
           color_as_factor = FALSE, background_color = "black")

spatPlot2D(gobject = visium_24,  show_image = TRUE, point_size = 4, 
           point_alpha = 10, cell_color = "noTLS_B", cell_color_gradient = custom_colors2,
           color_as_factor = FALSE)

#GENE EXPRESSION

HLA <- c("HHLA1", "HHLA2", "HHLA3", "HLA-A", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA",
         "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1",
         "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-E", "HLA-F", "HLA-G")

spatFeatPlot2D(gobject = visium_14, show_image = T, expression_values = "raw",
               feats = "HLA-DPA1", point_size = 2, cell_color_gradient = custom_colors3)
plotGiottoImage(gobject = all_samples, image_name = image_names)

