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
#use_python('C://python.exe', required = TRUE)
#use_miniconda(condaenv = 'C://r-miniconda/envs/r-reticulate/python.exe', required = TRUE)

use_miniconda(condaenv = "r-reticulate", required = TRUE)

#options(timeout = 6000000000) ### set this to avoid timeout error
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
#devtools::install_github("edsgard/trendsceek")
setwd("S:/VISIUM/basic_expression/")

#1.- DATA PREPARATION
results_folder = "S:/VISIUM/basic_expression/plots/"

#py_config()
my_python_path = 'C://python.exe'
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = TRUE,
                                  python_path = my_python_path)

#1.1.- open each of the samples
#ChIO SAMPLES
results_dir_40 <- "C://outs_p40/"
visium_40 = createGiottoVisiumObject(visium_dir = results_dir_40,
                                     expr_data = 'filter',
                                     png_name = 'tissue_lowres_image.png',
                                     gene_column_index = 2,
                                     instructions = instrs)


results_dir_50 <- "C://outs_p50/"
visium_50  = createGiottoVisiumObject(visium_dir = results_dir_50,
                                      expr_data = 'filter',
                                      png_name = 'tissue_lowres_image.png',
                                      gene_column_index = 2,
                                      instructions = instrs)


#1.2.- open and load the updated metadata onto each visium object
#these are the loupe-browser selected areas of interest. 
read.csv("C://metadata_40.csv", header = TRUE) -> metadata_40
addCellMetadata(gobject = visium_40, new_metadata = metadata_40, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_40
read.csv("C://metadata_50.csv", header = TRUE) -> metadata_50
addCellMetadata(gobject = visium_50, new_metadata = metadata_50, by_column = TRUE, column_cell_ID = "cell_ID") -> visium_50

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
#2.3.- actual filter
all_samples <- filterGiotto(gobject = all_samples,
                            expression_threshold = 1,
                            feat_det_in_min_cells = 5, #min num of spots (cell) that need to express a gene (feature)
                            min_det_feats_per_cell = 300, #min num of genes (features) that need to be detected in a spot (cell)
                            expression_values = c('normalized'),
                            verbose = T)
gc()

rm(list = c("metadata_40", "metadata_50"))
rm(list = c("visium_40", "visium_50"))
gc()

#trendSceek(all_samples)

#3.- SUBSET
# Extract logical indices for the subset, then add its corresponding metadata and split it from the merged samples
TLS_indices <- all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$TLS ==  1 #| all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B == 1
TLS_metadata <- all_samples@cell_metadata$cell$rna@metaDT[TLS_indices, ]
TLS_metadata$cell_ID -> TLS_spots
TLS_samples <- subsetGiotto(all_samples, cell_ids = TLS_spots)

mTLS_indices <- all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$mTLS ==  1 #| all_samples@cell_metadata[["cell"]][["rna"]]@metaDT$noTLS_B == 1
mTLS_metadata <- all_samples@cell_metadata$cell$rna@metaDT[mTLS_indices, ]
mTLS_metadata$cell_ID -> mTLS_spots
mTLS_samples <- subsetGiotto(all_samples, cell_ids = mTLS_spots)

saveGiotto(
  gobject = mTLS_samples ,
  foldername = "mTLSs_giotto",
  dir = "C:/",
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
  python_path = 'C://python.exe',
  init_gobject = TRUE,
  verbose = TRUE
) -> mTLStest

#pDataDT(all_samples) -> regulators
#4- DIFFERENTIAL EXPRESSION ANALYSIS
#4.1 SCRAN

table(numbercheck$response)
table(mTLS_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["response"]])
diff_expr_results1 <- findScranMarkers(gobject = mTLS_samples,
                                      expression_values = "normalized",
                                      cluster_column = "response", 
                                      group_1 = "CPR", 
                                      group_1_name = "CPR",
                                      group_2 = "NCPR", 
                                      group_2_name = "NCPR")
diff_expr_results1[[1]] -> CPR_POV
diff_expr_results1[[2]] -> NCPR_POV

#5.1. PLOTS
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

TLS_markers <- list(TLS_Van = c("CD3", "CD4", "CD8", "CD20", "MS4A1", "CD21", "CD23"), 
                    TLS_Cruz = c("MS4A1", "CD19", "CR2", "CXCL13"),
                    TLS_meylan = c("APOE", "C1QA", "C7", "CD52", "CD79A", "CXCL12", "DERL3", "FCRL5", "IGHA1", "IGHG1", 
                                   "IGHG2", "IGHG3", "IGHG4", "IGHGP", "IGHM", "IGKC", "IGLC1", "IGLC2", "IGLC3", 
                                   "IL7R", "JCHAIN", "LUM", "MZB1", "PIM2", "PTGDS", "PTLP", "SSR4", "TRBC2", "XBP1"),
                    TLS_mature = c("CD4", "PD1", "CXCR5", "CXCR13", "BCL6", "CD20", "MS4A1", "CD23", "PNAd"),
                    citokines = c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")
)

#heatmaps
plotMetaDataHeatmap(mTLS_samples, selected_feats = HLA, metadata_cols = c('response'), expression_values = "normalized",
                    x_text_size = 16, y_text_size = 16, strip_text_size = 20 )

ggsave(
  filename = ".tiff",
  plot = last_plot(),
  path = "S://",
  units = "in",
  width = 898 / 72,
  height = 1436 / 72, 
  dpi = 300 #increase to 300 if resolution needs to be corrected
)


violinPlot(mTLS_samples, feats = HLA, cluster_column = "response", strip_position = "top", 
           color_violin = "cluster", cluster_color_code = c("#0B645A", 'grey'), expression_values = "normalized",
           strip_text = 24, axis_text_x_size = 17, axis_text_y_size = 17)

ggsave(
  filename = ".tiff",
  plot = last_plot(),
  path = "C://",
  units = "in",
  width = 898 / 72,
  height = 2436 / 72, 
  dpi = 300 #increase to 300 if resolution needs to be corrected
)


#volcano plots
library(ggrepel)
volcano_data <- data.frame(gene = CPR_POV$feats, logFC = CPR_POV$summary.logFC, adjpval = -log10(CPR_POV$FDR), FDR = CPR_POV$FDR)
HLA_volcano_data <- subset(volcano_data, gene %in% HLA)
IG_volcano_data <- subset(volcano_data, gene %in% IG)

logFC_threshold <- 0.5
fdr_threshold <- 0.001

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
  filename = ".tiff",
  plot = last_plot(),
  path = "C:/",
  units = "in",
  width = 975 / 72,
  height = 1036 / 72, 
  dpi = 72 #increase to 300 if resolution needs to be corrected
)

rm(list = c("volcano_data", "NO_POV", "YES_POV", "diff_expr_results"))

##gene enrichment analysis
#open up DEG results (scran)
read.csv("C://mTLS_CPRpov.csv", header = 1, row.names = 1) -> de_genes

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
  filename = "GO_BP_CPR_mTLSs.tiff",
  plot = last_plot(),
  path = "C://",
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

write.csv(gsea_results_go, "C:/" )

#gsea_results_go@result -> gsea_df

#gsea plotting once its been calculated - revisiting results
#first open DEG genes, used as input for the GSEA
read.csv("C:/", header = 1, row.names = 1) -> de_genes

gene_ranking <- de_genes$summary.logFC
names(gene_ranking) <- de_genes$feats
gene_ranking <- sort(gene_ranking, decreasing = TRUE)
head(gene_ranking) #quickcheck

#then, load the GSEA results 
read.csv("C:/", header = TRUE) -> gsea

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
#after studying the GSEA results, can filter and order for choosing pathways to be plotted
selected_pathways_BP <- c("cell killing", "...",
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
  filename = ".tiff",
  plot = last_plot(),
  path = "C://",
  units = "in",
  width = 1138 / 72,
  height = 1236 / 72,   
  dpi = 300 #increase to 300 if resolution needs to be corrected
)


####TISSUE / SAMPLE PLOTTING
##############################
custom_colors3 = c('#FFF7F3', '#FDE0DD', '#FCC5C0', '#FA9FB5', '#F768A1', '#DD3497', '#AE017E', '#7A0177', '#49006A')
spots40 <- grep("^(40-)", all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]], value = TRUE)
#as samples were joined by shifting, this means the X axis was extended. Through the metadata, check the X range for the sample to ble plotted
sub40 <- subsetGiottoLocs(all_samples, x_max = 102761, x_min = 83427, return_gobject = TRUE) 
(names(all_samples@images)) -> image_names


#CELL-TYPE BINARY
spatPlot2D(gobject = sub40,  show_image = FALSE, point_size = 4,
                    cell_color = "mTLS", cell_color_code = c("grey", "#A81485"),#point_alpha = 10,
                    color_as_factor = FALSE, background_color = "white", title = "mTLSs",
                    point_shape = "border", point_border_col = "black", point_border_stroke = 0.3, 
                    legend_text = 12, axis_title = 29, axis_text = 19) +
  ggplot2::theme(plot.title = element_text(face = "bold", size = 45),
                 legend.text = element_text(size = 20))


#GENE EXPRESSION
HLA <- c("HHLA1", "HHLA2", "HHLA3", "HLA-A", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DOA",
         "HLA-DOB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1",
         "HLA-DQB2", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-E", "HLA-F", "HLA-G")

spatFeatPlot2D(gobject = visium_14, show_image = T, expression_values = "raw",
               feats = "HLA-DPA1", point_size = 2, cell_color_gradient = custom_colors3)
plotGiottoImage(gobject = all_samples, image_name = image_names)


