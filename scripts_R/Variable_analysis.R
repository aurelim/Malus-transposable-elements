library(scatterplot3d) 
library(FactoMineR)
library(factoextra)
library(sva)
library(gridExtra)
library(cowplot)
library(data.table)
library(tidyverse)
library(ggpubr)
library(FactoInvestigate)
library(cluster)
library(UpSetR)
library(ComplexUpset)
library(fuzzyjoin)
library(RColorBrewer)
options(digits = 4) # controls the number of significant digits to print
rm(list=ls())

source("scripts_R/fonctions.R")

table_name = fread("scripts_R/file_names.csv")
load("Antho_mappingTEdatabase_PCA/results/data_pca.rda")
table_name = filter(table_name, name %in% rownames(dt))
load("Antho_mappingTEdatabase_PCA/results/resPCA_TE.rda")
TE_class = fread("Antho_annotation_TEdatabase/Pastec_TE_sorter_Annotation_final_V3.csv", col.names = c("TE", "Order"))
  # mutate(TE = str_remove_all(string = TE, pattern = "\\|.*"))

TE_mapping_GDDH13 = fread("Aure_mappingTEdatabase_GDDH13/results/v2/g1.fasta.out.gff", col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"), sep = "\t") %>%
  mutate(seqid = str_replace_all(seqid, "chr0", "chr")) %>% 
  filter(seqid != "chr0") %>% 
  mutate(TE = str_remove_all(pattern = "Target [\"]Motif:", attributes)) %>% 
  mutate(TE = str_remove_all(pattern = "[\"]", TE)) %>% 
  separate(TE, into = c("TE", "V1", "V2"), sep = " ") %>% 
  # mutate(TE = str_remove_all(string = TE, pattern = "\\|.*")) %>%
  select(-c(V1, V2, attributes, source, type, score, phase)) %>% 
  select(GDDH13_seqid = seqid, GDDH13_start = start, GDDH13_end = end, TE) 

TE_mapping_Msiev = fread("Aure_mappingTEdatabase_Msieversii/results/Msieversii_haploid_v2.chr.fa.out.gff", col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"), sep = "\t") %>%
  mutate(TE = str_remove_all(pattern = "Target [\"]Motif:", attributes)) %>% 
  mutate(TE = str_remove_all(pattern = "[\"]", TE)) %>% 
  separate(TE, into = c("TE", "V1", "V2"), sep = " ") %>% 
  # mutate(TE = str_remove_all(string = TE, pattern = "\\|.*")) %>% 
  select(-c(V1, V2, attributes, source, type, score, phase)) %>% 
  select(Msiev_seqid = seqid, Msiev_start = start, Msiev_end = end, TE)

TE_mapping_Msyl = fread("Aure_mappingTEdatabase_Msylvestris/results/Msylvestris_haploid_v2.chr.fa.out.gff", col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"), sep = "\t") %>%
  mutate(TE = str_remove_all(pattern = "Target [\"]Motif:", attributes)) %>% 
  mutate(TE = str_remove_all(pattern = "[\"]", TE)) %>% 
  separate(TE, into = c("TE", "V1", "V2"), sep = " ") %>% 
  # mutate(TE = str_remove_all(string = TE, pattern = "\\|.*")) %>% 
  select(-c(V1, V2, attributes, source, type, score, phase)) %>% 
  select(Msyl_seqid = seqid, Msyl_start = start, Msyl_end = end, TE)


## contributions of TEs to PCA axes
contributions = resPCA_TE$var$contrib %>% as.data.frame() %>% rownames_to_column(var = "TE") %>% 
  filter(TE != "*") %>% 
  # mutate(TE = str_remove_all(string = TE, pattern = "\\|.*")) %>%
  # mutate(TE = str_replace_all(string = TE, pattern = "rnd_5_family_1611_DNA_ReconFamilySize22,", replacement = "rnd_5_family_1611_DNA_ReconFamilySize22")) %>% 
  mutate(sum = Dim.1 + Dim.2 + Dim.3 + Dim.6 + Dim.8) %>% 
  gather(key = "Axe", value = "contribution", Dim.1:sum) %>% 
  left_join(TE_class, by = "TE") %>% 
  mutate(TE_class = ifelse(Order == "LTR" | Order == "LINE" | Order == "SINE", "Class1", 
                           ifelse(Order == "Helitron" | Order == "Helitron|SINE" | Order == "TIR", "Class2", "noCat")))

# Plot TEs constribution distribution by axes, considering TE families
contributions %>%
  filter(Axe %in% c(paste0("Dim.", 1:3), "Dim.6", "Dim.8", "sum")) %>%
  filter(! is.na(Order)) %>% 
  filter(Order != "TRIM") %>% 
  filter(Order != "mixture") %>% 
  ggplot(aes(x = Axe, y = contribution, fill = Order))+
  scale_fill_manual(values = c("Helitron" = "#1B9E77",
                               "Helitron|SINE" = "blue",
                               "pararetrovirus" = "purple",
                               "LINE" = "#D95F02",
                               "LTR" = "#7570B3",
                               "MITE" = "#E7298A",
                               "noCat" = "#66A61E",
                               "SINE" = "#E6AB02",
                               "SSR" = "#A6761D",
                               "TIR" = "#666666"))+
  geom_boxplot()+
  scale_y_log10(limits = c(0.0000001,0.01))+
  ylab("TEs contribution (%)")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))

contributions %>%
  filter(Axe %in% c(paste0("Dim.", 1:3), "Dim.6", "Dim.8")) %>%
  group_by(Axe) %>%
  summarise(
    min = min(contribution),
    max = max(contribution),
    mean = mean(contribution),
    median = median(contribution),
    sd = sd(contribution),
    IQR = IQR(contribution)
  )

# Selection best contributing variable on 5 axes and sum of 5 axes

N = 200

data = contributions %>% 
  filter(Axe %in% c("Dim.1","Dim.2","Dim.3","Dim.6","Dim.8")) %>% 
  pivot_wider(names_from = Axe, values_from = contribution)

top_contributions <- bind_rows(
  lapply(c("Dim.1","Dim.2","Dim.3","Dim.6","Dim.8"), function(dim_col) {
    data %>%
      arrange(desc(.data[[dim_col]])) %>% # Utilisation correcte de la référence dynamique de colonne
      slice_head(n = N) %>%          # Prendre les 200 premières lignes
      mutate(Source_Dim = dim_col)     # Ajouter une colonne pour indiquer la dimension source
  })
)


# top_contributions = rbind(
#   # contributions %>% filter(Axe == "sum") %>% slice_max(contribution, n=N),
#   contributions %>% filter(Axe == "Dim.1") %>% slice_max(contribution, n=N),
#   contributions %>% filter(Axe == "Dim.2") %>% slice_max(contribution, n=N),
#   contributions %>% filter(Axe == "Dim.3") %>% slice_max(contribution, n=N),
#   contributions %>% filter(Axe == "Dim.6") %>% slice_max(contribution, n=N),
#   contributions %>% filter(Axe == "Dim.8") %>% slice_max(contribution, n=N)
# ) 

save(top_contributions, file = "200_top_contributors.rda")

top_contributions_on_genomes = top_contributions %>% 
  left_join(TE_mapping_GDDH13, by = "TE") %>% 
  left_join(TE_mapping_Msiev, by = "TE") %>% 
  left_join(TE_mapping_Msyl, by = "TE")

# ComplexUpset

d = top_contributions %>% 
  select(TE, Source_Dim, Order) %>%
  distinct(TE, Source_Dim, .keep_all = TRUE) %>% 
  mutate(presence = 1) %>% 
  pivot_wider(names_from = Source_Dim, values_from = presence, values_fill = 0)
Axes = colnames(d)[3:7]
d[Axes] = d[Axes] == 1

ComplexUpset::upset(d, Axes,
                    base_annotations = list(
                      'Intersection size'=intersection_size(
                        counts = TRUE,
                        mapping = aes(fill = Order),
                      ) + scale_fill_manual(values = c("Helitron" = "#1B9E77",
                                                       "LINE" = "#D95F02",
                                                       "LTR" = "#7570B3",
                                                       "MITE" = "#E7298A",
                                                       "noCat" = "#66A61E",
                                                       "SINE" = "#E6AB02",
                                                       "SSR" = "#A6761D",
                                                       "TIR" = "#666666"))
                    ), width_ratio = 0.1, sort_sets = FALSE)



