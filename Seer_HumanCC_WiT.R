library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(readxl)
library(stringr)
library(EnvStats)
library(ggrepel)
library(truncnorm)
library(scales)
library(naniar)
library(stats)
library(toolbox)
library(wrProteo)

#Load Seer data
#seer_CChuman_df = read_tsv("/Users/sorianos/Code/SomaLogic/Seer_Human_proteinGroups.txt")
seer_CChuman_df = read.delim("/Users/sorianos/Code/SomaLogic/Seer_Human_proteinGroups.txt")
#seer_CChuman_df3 = readMaxQuantFile("/Users/sorianos/Code/SomaLogic/","Seer_Human_proteinGroups.txt",separateAnnot = FALSE)

seer_CChuman_anno_df = read_tsv("/Users/sorianos/Code/SomaLogic/Seer_Human_annotation_NPs_separate.txt")
seer_CChuman_anno_df = seer_CChuman_anno_df %>% dplyr::select(Condition, Experiment)
colnames(seer_CChuman_anno_df)[colnames(seer_CChuman_anno_df) == "Experiment"] = "Sample_NP"

#Filter out contaminants in MaxQuant DDA data
seer_CChuman_df = seer_CChuman_df %>% filter(Potential.contaminant == "") %>%
                                      filter(Only.identified.by.site == "") %>%
                                      filter(Reverse == "")
                             
#Transform dataset into same format as mouse seer file (seer_KMC_df) - create features, sample, NP columns; prepare for WiT
seer_CChuman_df = seer_CChuman_df %>% pivot_longer(names_to = "Sample_NP",
                                                   values_to = "LFQ_Intensity",
                                                   cols = starts_with("LFQ")) %>%
                                      dplyr::select(c("Protein.IDs", "Protein.names","Gene.names","Sample_NP","LFQ_Intensity")) %>%
                                      mutate(Sample_NP = gsub("LFQ.intensity.","", Sample_NP)) %>%
                                      mutate(Log2_LFQ_Intensity = ifelse(LFQ_Intensity == 0, NA, log2(LFQ_Intensity))) %>% dplyr::select(-LFQ_Intensity) %>%
                                      filter(!is.na(Log2_LFQ_Intensity)) %>%
                                      mutate(Sample = substring(Sample_NP, 1, 3),
                                             NP = ifelse(grepl("NP1", Sample_NP), "NP1",
                                                          ifelse(grepl("NP2", Sample_NP), "NP2",
                                                                 ifelse(grepl("NP3", Sample_NP), "NP3",
                                                                        ifelse(grepl("NP4", Sample_NP), "NP4", "NP5")))))

#Merge in Case/Control annotation per Sample_NP
seer_CChuman_df = merge(seer_CChuman_df, seer_CChuman_anno_df, by = "Sample_NP")

#Write no sparsity filter version of data to saved table
write.table(seer_CChuman_df, "~/Code/SomaLogic/Seer_HumanCC_raw_noSF.tsv", sep = "\t", row.names = FALSE)

# #Find number of protein IDs (removing NP duplicates) in seer KMC data - BEFORE ANY FILTERS
# seer_CChuman_proteins = seer_CChuman_df %>% dplyr::select(c(`Protein IDs`,`Gene names`)) %>% unique() #%>% nrow()

# Valid Value Filtering
# Filter to those proteins present in at least 50% of at least one of the classes (Case, Control)
filt_seer_CChuman_df = seer_CChuman_df %>%
  group_by(Condition) %>%
  mutate(max_samples = n_distinct(Sample_NP)) %>% #Calculate number of samples per status
  ungroup() %>%
  group_by(Condition, `Protein.IDs`) %>%
  mutate(num_detected = n_distinct(Sample_NP), frac_detected = num_detected / max_samples) %>% #Calculate detection fraction of each prot ID, separated by status
  ungroup() %>%
  dplyr::select(c("Sample_NP", "Protein.IDs", "Condition", "frac_detected")) %>%
  filter(frac_detected >= 0.50) %>% #Keep prot IDs w/ frac_detection > 50%
  dplyr::select("Sample_NP", "Protein.IDs") %>%
  unique() %>%
  mutate(presence_tag = "sufficient") #Add column for tag indicating validity of prot ID (row)

#Add presence tag data and remove other rows/prot ID entries
filt_prot_class_presence = seer_CChuman_df %>%
  left_join(filt_seer_CChuman_df, by = c("Sample_NP","Protein.IDs")) %>%
  mutate(presence_tag = ifelse(is.na(presence_tag), FALSE,TRUE)) %>%
  filter(presence_tag) %>%
  mutate(Feature = paste(NP, `Protein.IDs`, sep = "|"))

#Select the final data table for the comparisons
seer_HumanCC_final_df = filt_prot_class_presence %>%
  dplyr::select(c("Sample_NP", "Condition", "Protein.IDs", "Protein.names", "Gene.names", "Feature", "Log2_LFQ_Intensity"))

#Flat 50% Sparsity Filter - not used
# num_samples = n_distinct(seer_CChuman_df$Sample)
#
# seer_HumanCC_final_df = seer_CChuman_df %>%
#   #filter(Log2_LFQ_Intensity > 0) %>%
#   mutate(Feature = paste(NP, Protein.IDs, sep = "_")) %>%
#   group_by(Feature) %>%
#   mutate(percent_detection = n_distinct(Sample)/num_samples) %>%
#   ungroup() %>%
#   filter(percent_detection >= 0.5)

#Find number of protein IDs (removing NP duplicates) in seer KMC data - AFTER sparsity filter, BEFORE WiT >1 filter
seer_CChuman_proteins = seer_HumanCC_final_df %>% dplyr::select(c(`Protein.IDs`,`Gene.names`)) %>% unique() #%>% nrow()

#Write 50% sparsity filter version of data to saved table
write.table(seer_HumanCC_final_df, "~/Code/SomaLogic/Seer_HumanCC_raw_50SF.tsv", sep = "\t", row.names = FALSE)

#Wilcox test
#Transform data to wide form for Wilcox test
seer_HumanCC_WiT_df = seer_HumanCC_final_df %>%
  group_by(Feature, Condition) %>%
  summarize(data = list(Log2_LFQ_Intensity), count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = c(data, count)) %>%
  filter(!is.na(count_Case)) %>% 
  filter(!is.na(count_Control))

#Wilcox test - store pvalue and difference
seer_HumanCC_WiT_res_df = seer_HumanCC_WiT_df %>%
  filter(!((count_Case == 1) & (count_Control == 1))) %>%
  rowwise() %>%
  mutate(wt_pval = wilcox.test(unlist(data_Case), unlist(data_Control))$p.value) %>%
  mutate(Diff = median(unlist(data_Case)) - median(unlist(data_Control))) %>% 
  ungroup() %>% 
  mutate(wt_pval_BH = p.adjust(wt_pval, method = "BH"))

seer_HumanCC_WiT_anno_df = seer_HumanCC_final_df %>% 
  dplyr::select(c("Feature", "Protein IDs", "Protein names", "Gene names")) %>% 
  distinct(Feature, .keep_all = TRUE)

seer_HumanCC_WiT_res_df = merge(seer_HumanCC_WiT_res_df, seer_HumanCC_WiT_anno_df, by = "Feature")

seer_HumanCC_WiT_res_df = seer_HumanCC_WiT_res_df %>% 
  mutate(p_val_sig = ifelse(wt_pval_BH < 0.05, "Significant", "Non"),
         proteins_sig = ifelse(wt_pval_BH < 0.05, `Gene names`, "")
  )

#Plot WiT pvalue distribution - before multiple testing correction
seer_HumanCC_WiT_pval_plot = seer_HumanCC_WiT_res_df %>%
  ggplot(aes(wt_pval)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_histogram(breaks = seq(0,1,0.05), fill = "#22c6c9", color = "black", alpha = 0.5) +
  labs(x = "p-value",
       y = "Aptamer Count",
       title = "Seer Proteograph: Mouse",
       subtitle = "50% Sparsity Filter") +
  theme(plot.subtitle = element_text(face = "italic"))
ggsave("seer_HumanCC_WiT_pval_plot.png")
seer_HumanCC_WiT_pval_plot

#Wilcoxon volcano Plot
seer_HumanCC_WiT_volc_plot = seer_HumanCC_WiT_res_df  %>%
  ggplot(aes(x = Diff, y = -log(wt_pval_BH, 10), group = p_val_sig)) +
  geom_point(aes(color = p_val_sig, text = `Gene names`), size = 2, alpha = 0.5) +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
  geom_vline(xintercept = -2, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = proteins_sig), max.overlaps = 20) +
  scale_color_manual(values = c("black", "red")) +
  labs(title = "Seer Proteograph: Human",
       x = "Log2(median_Case-median_Control)",
       y = "-log10(p-value)",
       subtitle = "50% Sparsity Filter") +
  theme_bw() +
  xlim(-5,5) +
  theme(legend.position = "none") +
  theme(plot.subtitle = element_text(face = "italic"))
ggsave("seer_HumanCC_WiT_volc_plot.png")
seer_HumanCC_WiT_volc_plot

seer_HumanCC_WiT_volc_plotly = ggplotly(seer_HumanCC_WiT_volc_plot, tooltip = "proteins_sig")
saveWidget(seer_HumanCC_WiT_volc_plotly, "seer_HumanCC_WiT_volc_plotly.html", selfcontained = F, libdir = "lib")

#Remove list-type columns before writing WiT result table
w_seer_HumanCC_WiT_res_df = seer_HumanCC_WiT_res_df %>% dplyr::select(-c("data_Case","data_Control"))
write.table(w_seer_HumanCC_WiT_res_df, "~/Code/SomaLogic/Seer_HumanCC_WiT.tsv", sep = "\t", row.names = FALSE)


