##########################################
## 1.- ecology on deconvolution results ##
##########################################

read.csv("S:/VISIUM/deconvolution/n15/n15_SpaCET_response_resultsv2.csv", header = TRUE) -> clin_nonmalignant
clin_nonmalignant[clin_nonmalignant$response == "NCPR",] -> NCPR_nonmalignant
NCPR_nonmalignant[!NCPR_nonmalignant$list_ID %in% c(4, 46, 72),] -> NCPR_nonmalignant

nonmalignant <- NCPR_nonmalignant[,c(-2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13 )]

library(vegan)
nonmalignant2 <- nonmalignant[,c(-1)]
rownames(nonmalignant2) <- nonmalignant[,1]
shannon_diversity <- diversity(nonmalignant2, index = "shannon")
richness <- specnumber(nonmalignant2)
evenness <- shannon_diversity / log(richness)

#14-AAACAAGTATCTCCCA-1 has a richness of 26, lets check it
sum(nonmalignant2["14-AAACAAGTATCTCCCA-1", ] > 0) #indeed



diversity_df <- data.frame(cell_ID = names(shannon_diversity),
                           shannon_diversity = shannon_diversity,
                           richness = richness,
                           evenness = evenness)


read.csv("S:/VISIUM/TFM/15 samples/final_metadata_NCPRs.csv")[,-1] -> metadata_NCPR

merge(metadata_NCPR, diversity_df, by = "cell_ID") -> Experimental

write.csv(diversity_df, "S:/VISIUM/TFM/15 samples/ecology/deconvolution_ecology_NCPRs.csv")

library(tidyverse)

Experimental %>% relocate(area, .after = list_ID) -> Experimental
Experimental <- Experimental[,c(-2)]

Experimental[Experimental$area == "mTLSs",] -> mTLSs
Experimental[Experimental$area == "Low B cell areas",] -> lowBcellareas

hist(lowBcellareas$richness)



long_data5 <- Experimental %>%
  pivot_longer(cols = "evenness", 
               values_to = "Expression") %>%
  mutate(evenness = factor(evenness))

long_data5$group <- interaction(long_data5$Tumor_stroma_interface, long_data5$progression, sep = "_")
long_data5$progression <- factor(long_data5$progression, levels = c("YES", "NO"))

long_data6 <- Experimental %>%
  pivot_longer(cols = "shannon_diversity", 
               values_to = "Expression") %>%
  mutate(shannon_diversity = factor(shannon_diversity))

long_data6$group <- interaction(long_data6$Tumor_stroma_interface, long_data6$progression, sep = "_")
long_data6$progression <- factor(long_data6$progression, levels = c("YES", "NO"))

long_data7 <- Experimental %>%
  pivot_longer(cols = "richness", 
               values_to = "Expression") %>%
  mutate(richness = factor(richness))

long_data7$group <- interaction(long_data7$Tumor_stroma_interface, long_data7$progression, sep = "_")
long_data7$progression <- factor(long_data7$progression, levels = c("YES", "NO"))


ggplot(long_data7, aes(
  x = factor(`Tumor_stroma_interface`, levels = c("Tumor", "Interface", "Stroma")), 
  y = Expression,
  fill = factor(progression, levels = c("YES","NO"))
)) +
  geom_jitter(shape = 16, size = 1.5, color = "grey45", width = 0.4, alpha = 0.6, aes(group = interaction(Tumor_stroma_interface, progression))) + 
  geom_violin(aes(group = interaction(Tumor_stroma_interface, progression)), trim = TRUE, linewidth = 1, position = position_dodge(width = 0.7)) +
  stat_summary(aes(group = interaction(Tumor_stroma_interface, progression)), fun = median, geom = "point", size = 2, color = "black", position = position_dodge(width = 0.7)) +
  # scale_y_log10() +
  labs(x = "", y = "Richness", title = "", fill = "") +
  scale_fill_manual(values = c("#C0C0FF", "#FFE0E0")) + # Custom fill colours
  scale_colour_manual(values = c("#000080", "#A00000")) + #custom line colours # Customize progression colors
  theme_minimal() + 
  theme(
    axis.text.y = element_text(size = 22),
    axis.text.x = element_text(size = 24, face = "bold", angle = 30, hjust = 1),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
    legend.text = element_text(size = 22),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 20),
    legend.position = "bottom",
    strip.text = element_text(size = 30, face = "bold"),
    panel.spacing = unit(1.5, "lines"),
    legend.key.size = unit(10, "mm")
  ) +
  guides(colour = "none")


ggsave(                      
  filename = "richness_deconvolution.tiff",
  plot = last_plot(),
  path = "S:/VISIUM/TFM/15 samples/ecology/",
  units = "in",
  width = 1000 / 72, 
  height = 536 / 72,   
  dpi = 300 #72 or 300
)

#########################################
###### 2.- ecology on L-R results #######
#########################################

read.csv("S:/VISIUM/TFM/15 samples/LRs/aggegated_LRs_expr.csv", header = TRUE, row.names = 1) -> aggregated_LR
read.csv("S:/VISIUM/TFM/15 samples/LRs/LR_filters_agg.csv", header = TRUE, row.names = 1) -> LR_filters_agg
read.csv("S:/VISIUM/TFM/15 samples/LRs/per_sample_agg.csv", header = TRUE, row.names = 1) -> per_sample_agg

#with the aggregated_LR dataframe

#1.-Compute per-spot diversity (Shannon), richness and evenness using product_expr as abundance.
library(data.table)
DT <- as.data.table(aggregated_LR)
DT[, product_expr := as.numeric(product_expr)]

# Compute per-spot summaries: Shannon (natural log), richness, evenness
# We treat product_expr as abundance; convert to relative frequency per spot
# and compute Shannon: -sum(p * log(p)), ignoring p==0.
#first, excludes product expression below 0. Groups are filtered by spot and sample
#per spot it calculates the sum of the product expression and the richness (number of LRs)
spot_div <- DT[product_expr > 0,                    # ignore zeros for p*log(p)
               .(sum_prod = sum(product_expr),
                 richness = uniqueN(LR_pair, na.rm = TRUE)), # temp: counts of LR pairs with >0 per spot
               by = .(sample_id, spot_id)]             # group for each spot

#second step, conceptual one: the product expr per spot.
#proportions are calculated as product expression in spot/sum of all product expr values in the spot
# But we need Shannon, so join back and compute p
# Alternative single-pass:
spot_div <- DT[, {
  v <- product_expr
  # if all zeros -> Shannon = 0, richness = 0
  if (all(is.na(v)) || sum(v, na.rm=TRUE) == 0) {
    list(shannon = 0, richness = 0, sum_prod = 0) #there are no LR in a spot where the expression product is 0, neither spots with no LRs in this table
  } else {
    vpos <- v[v > 0]
    p <- vpos / sum(vpos) #proportions are calculated as product expression in spot/sum of all product expr values in the spot
    sh <- -sum(p * log(p))                   # natural log; log2 if you want bits
    list(shannon = sh, richness = length(vpos), sum_prod = sum(v, na.rm=TRUE))
  }
}, by = .(sample_id, spot_id, progression, Tumor_stroma_interface)]

#in a spot, you take an LR's expression product
#that is relative to that spot total of LRs expression products sum, being the relative frequency
#check it is correct
# mtls spot 49-30x68 has 3 LRs with product expression of 77.79979
# LR1 = APP_CD74; PRexp = 19.36517; p = 0.24891; p * log(p) = -0.3461501
# LR2 = CD36_COL1A1; PRexp = 26.13631; p = 0.33594; p * log(p) = -0.366451
# LR3 = CD36_COL1A2; PRexp = 32.29831; p = 0.41515; p * log(p) = -0.3649647
#shannon diversity is 1.077566 (my calculation) // 1.077567 (looped)
#IT IS CORRECT (i did not take all decimals)
#how many spots have at least one LR? 31901 out of 42797 for 12 samples
#in the 15 samples - NCPR group, it is 20536 spots out of 28701
unique(DT$spot_id) -> spots_withLRs 
#the results of these must be joined to the total, so we can see the spots with 0 LRs per Tumor_stroma_interface


# Evenness = shannon / log(richness) (if richness>1)
spot_div[, evenness := ifelse(richness > 1, shannon / log(richness), NA)]

write.csv(spot_div, "S:/VISIUM/TFM/15 samples/ecology/LRs_diversity.csv")
# Inspect
head(spot_div)

library(tidyverse) 
#plot
plot_df <- spot_div %>%
  mutate(
    progression = factor(progression, levels = c("YES", "NO")),
    Tumor_stroma_interface = factor(Tumor_stroma_interface, levels = c("Tumor", "Stroma", "Interface"))
  )

ggplot(plot_df, aes(
  x = Tumor_stroma_interface,
  y = shannon, #evenness, richness, shannon
  fill = progression
)) +
  geom_jitter(
    shape = 16, size = 1.5, color = "grey45",
    width = 0.4, alpha = 0.6,
    aes(group = interaction(Tumor_stroma_interface, progression))
  ) +
  geom_violin(
    aes(group = interaction(Tumor_stroma_interface, progression)),
    trim = TRUE, linewidth = 1,
    position = position_dodge(width = 0.7)
  ) +
  stat_summary(
    aes(group = interaction(Tumor_stroma_interface, progression)),
    fun = median, geom = "point", size = 2, color = "black",
    position = position_dodge(width = 0.7)
  ) +
  labs(
    x = "", y = "Shannon diversity",
    title = "", fill = ""
  ) +
  scale_fill_manual(values = c("#C0C0FF", "#FFE0E0")) +
  scale_colour_manual(values = c("#000080", "#A00000")) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 22),
    axis.text.x = element_text(size = 24, face = "bold", angle = 30, hjust = 1),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
    legend.text = element_text(size = 22),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 20),
    legend.position = "bottom",
    strip.text = element_text(size = 30, face = "bold"),
    panel.spacing = unit(1.5, "lines"),
    legend.key.size = unit(10, "mm")
  ) +
  guides(colour = "none")

ggsave(
  filename = "LRs_shannon_div.tiff",
  plot = last_plot(),
  path = "S:/VISIUM/TFM/15 samples/ecology/",
  units = "in",
  width = 1000 / 72,
  height = 536 / 72,
  dpi = 300
)


###### ecological statistics ######
read.csv("S:/VISIUM/TFM/15 samples/final_metadata_NCPRs.csv", header = TRUE, row.names = 1) -> visium_metadata #1:22
read.csv("S:/VISIUM/TFM/15 samples/ecology/deconvolution_ecology_NCPRs.csv", header = TRUE, row.names = 1) -> deconv_ecology
colnames(deconv_ecology) <- c("cell_ID", "shannon_deconv", "richness_deconv", "evenness_decnov")
read.csv("S:/VISIUM/TFM/15 samples/ecology/LRs_diversity.csv", header = TRUE, row.names = 1)[,-c(1, 3, 4)] -> LRs_ecology
colnames(LRs_ecology) <- c("spot_id",  "shannon_LRs", "richness_LRs", "sum_prod", "evenness_LRs")
 
visium_metadata$spot_id <- paste0(substr(visium_metadata$cell_ID, 1, 3), visium_metadata$spot_id) #give sample prefix to spot_id
#create a dataframe that matches the spacet and 10x barcode IDs
visium_metadata[, c(1,2)] -> IDs_key
merge(IDs_key, LRs_ecology, by = "spot_id") -> LRs_ecology
merge(deconv_ecology, LRs_ecology, by = "cell_ID", all.x = TRUE) -> ecological_data

merge(visium_metadata, ecological_data, by = "cell_ID") -> complete_data
complete_data <- complete_data[, -c(23, 27)]

#write.csv(complete_data, "S:/VISIUM/TFM/15 samples/ecology/ecology_scores.csv" )

#calculate statistics per area of interest
complete_data[complete_data$Tumor_stroma_interface == "Tumor",] -> diversity #4286
complete_data[complete_data$Tumor_stroma_interface == "Stroma",] -> diversity #21343
complete_data[complete_data$Tumor_stroma_interface == "Interface",] -> diversity #3072



cal3 = colnames(diversity)[23:length(colnames(diversity))]
RE_3 = as.data.frame(cal3)
colnames(RE_3) = c("Variable")

# Genera nueva columna titulada P_VALUE
RE_3$P_VALUE = c(NA)
# Genera nuevas columnas para calcular otras variables
RE_3$P_VALUE_ADJ = c(NA)
RE_3$Median_PD = c(NA)
RE_3$Median_NPD = c(NA)
RE_3$Mean_PD = c(NA)
RE_3$Mean_NPD = c(NA)
RE_3$N_PD = c(NA)
RE_3$N_NPD = c(NA)
RE_3$P25_PD = c(NA)
RE_3$P75_PD = c(NA)
RE_3$P25_NPD = c(NA)
RE_3$P75_NPD = c(NA)

#table(Experimental$response)
#write.csv(df_new,"S:/VISIUM/deconvolution/analysis/RESULTS/no_malignant/SpaCET_nomalignant.csv", row.names = TRUE, col.names=TRUE)
#write.csv(final,"S:/VISIUM/deconvolution/analysis/RESULTS/no_malignant/SpaCET_nomalignant_RESPONSE.csv", row.names = TRUE, col.names=TRUE)

# Ejecuta analisis U Mann Whitney
# Cuando haya un número entero, desde la columna 6 (donde empiezan las variables) hasta el número último de columnas
for (l in 23:ncol(diversity)) {
  x = wilcox.test(diversity[diversity$progression == "YES", l],
                  diversity[diversity$progression == "NO", l])
  
  pv = x$p.value
  RE_3$P_VALUE[l-22] = pv
  RE_3$N_NPD[l-22] = length(na.omit(diversity[diversity$progression == "NO", l]))
  RE_3$N_PD[l-22] = length(na.omit(diversity[diversity$progression == "YES", l]))
  RE_3$Median_NPD[l-22] = median(diversity[diversity$progression == "NO", l], na.rm = TRUE)
  RE_3$Median_PD[l-22] = median(diversity[diversity$progression == "YES", l], na.rm = TRUE)
  RE_3$Mean_NPD[l-22] = mean(diversity[diversity$progression == "NO", l], na.rm = TRUE)
  RE_3$Mean_PD[l-22] = mean(diversity[diversity$progression == "YES", l], na.rm = TRUE)
  RE_3$P25_NPD[l-22] = quantile(diversity[diversity$progression == "NO", l], 0.25, na.rm = TRUE)
  RE_3$P75_NPD[l-22] = quantile(diversity[diversity$progression == "NO", l], 0.75, na.rm = TRUE)
  RE_3$P25_PD[l-22] = quantile(diversity[diversity$progression == "YES", l], 0.25, na.rm = TRUE)
  RE_3$P75_PD[l-22] = quantile(diversity[diversity$progression == "YES", l], 0.75, na.rm = TRUE)
}


# Evitar que los decimales salgan con notación científica
RE_3$Median_NPD = format(RE_3$Median_NPD, scientific = FALSE)
RE_3$Median_PD = format(RE_3$Median_PD, scientific = FALSE)
RE_3$Mean_NPD = format(RE_3$Mean_NPD, scientific = FALSE)
RE_3$Mean_PD = format(RE_3$Mean_PD, scientific = FALSE)
RE_3$P25_PD = format(RE_3$P25_PD, scientific = FALSE)
RE_3$P25_NPD = format(RE_3$P25_NPD, scientific = FALSE)
RE_3$P75_PD = format(RE_3$P75_PD, scientific = FALSE)
RE_3$P75_NPD = format(RE_3$P75_NPD, scientific = FALSE)

padj = p.adjust(RE_3$P_VALUE, method= "fdr", n = length(RE_3$P_VALUE))
RE_3$P_VALUE_ADJ = padj

# Guarda los resultados
write.csv(RE_3, "S:/VISIUM/TFM/15 samples/ecology/diversity_Tumor_PDvsNPD.csv", row.names = TRUE, col.names=TRUE)



### bonus, some top L-Rs per region
################################################################################

top_15_per_group <- LR_filters_agg %>%
  group_by(Tumor_stroma_interface, progression) %>%
  slice_max(order_by = median_product_expr, n = 15) %>%
  ungroup()

write.csv(top_15_per_group, "S:/VISIUM/TFM/Res_Sofia/my_launch/nocells/top15_LRs.csv")

mean_product_by_progression <- all_LR_keep_ud_exp %>%
  group_by(LR_pair, Tumor_stroma_interface, progression) %>%
  summarise(
    mean_product_expr = mean(product_expr, na.rm = TRUE),
    .groups = "drop"
  )

delta_LR_expr <- mean_product_by_progression %>%
  pivot_wider(
    names_from = progression,
    values_from = mean_product_expr,
    values_fill = 0
  ) %>%
  mutate(
    delta_product_expr = NO - YES
  )

deltatop_15 <- delta_LR_expr %>%
  group_by(LR_pair, Tumor_stroma_interface) %>%
  slice_max(order_by = delta_product_expr, n = 15) %>%
  ungroup()
write.csv(deltatop_15, "S:/VISIUM/TFM/Res_Sofia/my_launch/nocells/deltatop15_LRs.csv")

################################################################################