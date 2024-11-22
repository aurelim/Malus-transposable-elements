library(tidyverse)
library(data.table)

options(digits = 4) # controls the number of significant digits to print
rm(list=ls())

source("scripts_R/fonctions.R")

file_names <- dir("Antho_mappingTEdatabase/results/", pattern = ".txt")
table_name = fread("scripts_R/file_names.csv")
file_names = paste0(table_name$name, "_list_sorting.txt")
TE_names = fread(paste0("Antho_mappingTEdatabase/results/",file_names[1]))


# step 1 : As usual, normalisation sur tous les individus

list <- vector("list", length(file_names))

for (i in 1:length(file_names)) {
  file <- read.table(paste0("Antho_mappingTEdatabase/results/",file_names[i]), header = FALSE, sep = "\t")
  #file <- subset(file, V2 >1000)
  total <- sum(file$V3) #Nb de reads mappés total sur la base de données
  normalized <- (file$V3 * table_name$reads[i]) / (file$V2 * 0.001 * total  * 0.000001)
  list[[i]] <- normalized
}
names(list) <- str_remove_all(file_names, pattern = "_list_sorting.txt")
df <- as.data.frame(list)


# step 2 : CombatSeq https://github.com/zhangyuqing/ComBat-seq?tab=readme-ov-file

data.active = df
data.active[is.na(data.active)] <- 0
batch = as.data.frame(colnames(data.active)) %>% full_join(table_name, by = c("colnames(data.active)" = "name")) %>% 
  mutate(code = ifelse(who == "Xilong", 1, 2)) %>% select(code)

adjusted <- ComBat_seq(data.active, batch=batch$code, group=NULL)

save(adjusted, file = "Antho_mappingTEdatabase_PCA/results/adjusted.rda")
load("Antho_mappingTEdatabase_PCA/results/adjusted.rda")



data_PCA = t(adjusted)
colnames(data_PCA) = TE_names$V1

# Remove some individuals (OTH and Baccata)

list_OTH = (filter(table_name, specie == "OTH"))$name
list_baccata = (filter(table_name, specie == "baccata"))$name
data_PCA = data_PCA[! rownames(data_PCA) %in% list_OTH, ]
data_PCA = data_PCA[! rownames(data_PCA) %in% list_baccata, ]
data_PCA = data_PCA[! rownames(data_PCA) %in% c("Dom_09", "Dom_22", "Syl_01", "Syl_07"), ]
data_PCA = data_PCA[! rownames(data_PCA) %in% c("Dom_15", "Dom_29", "Dom_35", "Dom_34", "Dom_32", "Syl_10", "D12"),]

table_name = filter(table_name, name %in% rownames(data_PCA))

unique(table_name$specie)

data_PCA = as.data.frame(data_PCA)

dt = left_join(rownames_to_column(data_PCA[,1:158332]), table_name[,c(1,3,4)], by=c("rowname" = "name")) %>% column_to_rownames(var = "rowname") %>% 
  filter(specie != "baccata")

save(dt, file = "Antho_mappingTEdatabase_PCA/results/data_pca.rda")
