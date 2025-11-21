library(Giotto)
library(tidyverse)
library(SpaCET)

##load the "all_samples" normalized object and then split by sample
path_all_samples <- "S:/VISIUM BMS/MARTA_BMS/allsamplesBMS_giotto"
all_samples <- loadGiotto(path_all_samples,
                          load_params = list(),
                          reconnect_giottoImage = TRUE,
                          python_path = NULL,
                          init_gobject = TRUE,
                          verbose = TRUE
)

#check the samples available
#add Marta's metadata
read.csv("S:/VISIUM BMS/MARTA_BMS/metadata_filtered.csv", header = TRUE) -> metadata_MMA_BMS
addCellMetadata(gobject = all_samples, new_metadata = metadata_MMA_BMS, by_column = TRUE, column_cell_ID = "cell_ID") -> all_samples
pDataDT(all_samples) -> metadata
print(unique(metadata$list_ID))


#split giotto objects per sample
samples <- c(12, 21, 23, 27, 30, 37, 45)
for (id in samples) {
  # extract matching spots
  spots <- grep(paste0("^(", id, "-)"),
                all_samples@cell_metadata[["cell"]][["rna"]]@metaDT[["cell_ID"]],
                value = TRUE)
  
  # subset Giotto object
  sub <- subsetGiotto(all_samples, cell_ids = spots)
  
  # assign to environment
  assign(paste0("spots", id), spots)
  assign(paste0("sub", id), sub)
}

rm(sub)


######## SPACET
#1. GENERATE SpaCET OBJECT
# it does not matter which visium object you load as you will only use the SpaCET object architecture (substituting its data with your own), regardless of the sample
#you use it a an skeleton

visiumPath <- "S:/VISIUM/visium_results/outs_p40/"
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)

#2. LOAD YOUR DATA TO THE SpaCET OBJECT

# each sample will be done individually; change "subXX" as needed

sub30@expression[["cell"]][["rna"]][["normalized"]]@exprMat -> SpaCET_obj@input[["counts"]]
sub30@expression[["cell"]][["rna"]][["normalized"]]@exprMat@Dimnames -> SpaCET_obj@input[["counts"]]@Dimnames
sub30@expression[["cell"]][["rna"]][["normalized"]]@exprMat@x -> SpaCET_obj@input[["counts"]]@x

#3.- SpaCET DECONVOLUTION: WATCH OUT FOR TUMOUR HISTOLOGY
#check the sample in order to execute the matching deconvolution

### Squamous · LUSC 
# patients: 2 · 14 · 23 · 24 · 37 · 12 · 35 · 45 · 21 · 33
#for this analysis: 12 · 21 · 23 · 37 · 45

### Adeno · LUAD
# patients: 20 · 11 · 13 · 15 · 20 · 22 · 17 · 27 · 32 · 30 · 28 · 43 · 49
#for this analysis: 27 · 30


# deconvolve ST data
#execute one block, LUAD or LUSC

spacet_LUAD <- SpaCET.deconvolution(SpaCET_obj = SpaCET_obj, cancerType="LUAD", coreNo=8)
spacet_LUAD@results$deconvolution$propMat -> LUAD_deconv
t(LUAD_deconv) -> LUAD_deconv_proc

spacet_LUSC <- SpaCET.deconvolution(SpaCET_obj = SpaCET_obj, cancerType="LUSC", coreNo=8) 
spacet_LUSC@results$deconvolution$propMat -> LUSC_deconv
t(LUSC_deconv) -> LUSC_deconv_proc

# change the saved file name so they do not overwrite
write.csv(LUAD_deconv_proc, "S:/VISIUM BMS/MARTA_BMS/deconvolution_raw/p30_LUAD_deconv.csv", col.names = TRUE, row.names = TRUE)
write.csv(LUSC_deconv_proc, "S:/VISIUM BMS/MARTA_BMS/deconvolution_raw/p45_LUSC_deconv.csv", col.names = TRUE, row.names = TRUE)

table(metadata_MMA_BMS$list_ID)

################################
#PROCESSING
#once all samples have been deconvolved with spacet:

library(tidyverse)

path <- "S:/VISIUM BMS/MARTA_BMS/deconvolution_raw/" #wherever you have the spacet csv results
file_list <- list.files(path, pattern = "*.csv", full.names = TRUE)
final_df <- file_list %>%
  lapply(read.csv) %>%
  bind_rows()

colnames(final_df)[1] <- "cell_ID"
write.csv(final_df, "S:/VISIUM BMS/MARTA_BMS/deconvolution_raw/merged_deconvolution.csv" , col.names = TRUE, row.names = FALSE)




