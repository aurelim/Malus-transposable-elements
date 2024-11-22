library(scatterplot3d) 
library(FactoMineR)
library(factoextra)
# # install.packages("devtools")
# devtools::install_github("zhangyuqing/sva-devel", force = TRUE)
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

setwd("/home/aurelie/Dropbox/Post-Doc_ECLECTIC/2_Projet_TE_malus/")
source("scripts_R/fonctions.R")
table_name = fread("Antho_mappingTEdatabase_PCA/scripts/file_names.csv")
load("Antho_mappingTEdatabase_PCA/results/data_pca.rda")
table_name = filter(table_name, name %in% rownames(dt))

resPCA_TE <- PCA(dt, scale.unit = T, graph = F, ncp = 30, quali.sup = c(158334,158333))
save(resPCA_TE, file = "Antho_mappingTEdatabase_PCA/results/resPCA_TE.rda")
load("Antho_mappingTEdatabase_PCA/results/resPCA_TE.rda")

fviz_pca_ind(resPCA_TE, col.ind = dt$specie)+
  geom_point(aes(color = factor(dt$specie),
                 shape = factor(dt$specie2),
  ), 
  size = 2)+
  guides(shape = guide_legend(title = "species"),
         colour = guide_legend(title = "population"))+
  scale_shape_manual(values=c(rep(15,6), 16, 17, 18, 3, 4, 8, 11, 14,12, 15, 15))

fviz_eig(resPCA_TE, addlabels = TRUE, choice = "variance", ncp = 30 )

get_eig(resPCA_TE) %>% as.data.frame() %>% rownames_to_column(var = "Dim") %>%
  mutate(Dim = str_remove_all(Dim, "Dim.")) %>% 
  mutate(Dim = as.numeric(Dim)) %>% 
  ggplot(aes(x = Dim, y = cumulative.variance.percent))+
  geom_bar(stat = "identity")+
  xlim(c(0,30))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90))


ind_coord <- resPCA_TE[["ind"]][["coord"]]

colors <- c("red", "blue", "green", "purple", "orange", "gold")
colors <- colors[as.numeric(as.factor(table_name$specie))]

# shapes = c(16, 17)
# shapes <- shapes[as.numeric(as.factor(table_name$who))]

s3d_all <- scatterplot3d(ind_coord[,1:3], pch = 16, color=colors, grid=TRUE, box=FALSE)
legend("topright", pch=19, col=c("red", "blue", "green", "purple", "orange", "gold"), 
       legend=c("domestica", "domestica_r" ,"orientalis","sieversii", "sieviersii_k", "sylvestris"),
       cex= 1)
title(main="PCA TE")

A <- fviz_pca_ind(resPCA_TE,
                  geom.ind = "point", # show points only (but not "text")
                  pointshape = 21,
                  axes = c(1,2),
                  pointsize = 2,
                  fill.ind = table_name$specie,
                  col.ind = "black",
                  palette = c("red", "blue", "black", "green", "purple", "pink", "gold", "brown"),
                  addEllipses = FALSE, # Concentration ellipses
                  legend.title = "Groups",
                  mean.point = FALSE,
                  label = "ind"
)
B <- fviz_pca_ind(resPCA_TE,
                  geom.ind = "point", # show points only (but not "text")
                  pointshape = 21,
                  axes = c(1,3),
                  pointsize = 2,
                  fill.ind = table_name$specie,
                  col.ind = "black",
                  palette = c("red", "blue", "black", "green", "purple", "pink", "gold", "brown"),
                  addEllipses = FALSE, # Concentration ellipses
                  legend.title = "Groups",
                  mean.point = FALSE
)
C <- fviz_pca_ind(resPCA_TE,
                  geom.ind = "point", # show points only (but not "text")
                  pointshape = 21,
                  axes = c(2,3),
                  pointsize = 2,
                  fill.ind = table_name$specie,
                  col.ind = "black",
                  palette = c("red", "blue", "black", "green", "purple", "pink", "gold", "brown"),
                  addEllipses = FALSE, # Concentration ellipses
                  legend.title = "Groups",
                  mean.point = FALSE
)
D <- fviz_pca_ind(resPCA_TE,
                  geom.ind = "point", # show points only (but not "text")
                  pointshape = 21,
                  grid = TRUE,
                  axes = c(1,4),
                  pointsize = 2,
                  fill.ind = table_name$specie,
                  col.ind = "black",
                  palette = c("red", "blue", "black", "green", "purple", "pink", "gold", "brown"),
                  addEllipses = FALSE, # Concentration ellipses
                  legend.title = "Groups",
                  mean.point = FALSE
)


p = plot_grid(A + theme(legend.position = "none"), 
              B + theme(legend.position = "none"), 
              C + theme(legend.position = "none"), 
              D + theme(legend.position = "none"), 
              labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2)

legend <- get_legend(
  # create some space to the left of the legend
  A + theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(p, legend, rel_widths = c(2,.5))
