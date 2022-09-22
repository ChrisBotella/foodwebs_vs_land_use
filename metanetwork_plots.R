rm(list=ls())
#This script aims at visualising and comparing networks for Botella paper
#using metanetwork package

setwd("~/Bureau/article_chris_b/github/")
# mainDir ="~/Desktop/collab_LECA/article_chris_b/github/" 
# setwd(mainDir)

### metanetwork Installation
# library(devtools)
# Download package repository from
# https://gitlab.com/marcohlmann/metanetwork
# then install from source using devtools::install
require(metanetwork)
require(igraph)
require(ggplot2)
require(network)
require(intergraph)
library(stringr)
library(RColorBrewer)
library(Matrix.utils)
library(raster)
library(sp)
library(ggraph)
library(tidyverse)
library(scatterpie)
library(viridisLite)

# loading data
load("data/preprocessed_data.Rdata")

#building the igraph metaweb
metaweb = graph_from_adjacency_matrix(t(metaweb_Adj), mode = "directed")
cellsVsSpecies = cellsVsSpecies[rownames(cells),]
cellsVsSpecies_l = cellsVsSpecies

# build abundance table per landUse or per region or per intens
covariate = "intens"
rownames(cellsVsSpecies_l) = cells[,covariate]
abTable_loc = t(sapply(by(as.matrix(cellsVsSpecies_l),rownames(cellsVsSpecies_l),colMeans),identity))
#site names should not be a numeric
rownames(abTable_loc) = paste0(substr(covariate,1,2),rownames(abTable_loc))

#building trophicTable !!CAUTION NA in group at line  513!!!!!
spTab_loc = spTab[,c("name","SBMgroup","Class")]
rownames(spTab_loc) = spTab_loc$name
colnames(spTab_loc) = c("species","SBMgroup","Class")
spTab_loc[513,"SBMgroup"] = 16 #assigning NA species to arbitrary group

#adding diets to trophicTable
diets = V(metaweb)$name[ which(!(V(metaweb)$name %in% rownames(spTab_loc)))]
spTab_loc = rbind(spTab_loc,cbind(species = diets,SBMgroup = diets,Class = "Basal resources"))
rownames(spTab_loc) = spTab_loc$species
#adding diets to abTable (uniform abundance, equals to )
value_diets = min(abTable_loc[which(abTable_loc>0)])
df = data.frame(matrix(rep(value_diets,length(diets)),nrow = 1)) 
df = rbind(df, df[rep(1, nrow(abTable_loc)-1),])
colnames(df) = diets
rownames(df) = rownames(abTable_loc)
abTable_loc = cbind(abTable_loc,df)

#load saved meta_verte object (with layout)
load("foodwebs_vs_land_use/data/meta_verte.Rdata")

# # build metanetwork
# meta_verte = build_metanet(metaweb = metaweb,abTable = abTable_loc,
#                              trophicTable = spTab_loc)
# meta_verte = append_agg_nets(meta_verte)
# meta_verte = compute_TL(meta_verte)
# 
# metanetwork::print(meta_verte)

#METANETWORK ANALYSIS#
######################


#at a species level
##################

#attaching layout
beta = 0.00005
# meta_verte = attach_layout(g = meta_verte$metaweb,
#                          metanetwork = meta_verte,beta = beta)
# save(meta_verte,file="foodwebs_vs_land_use/data/meta_verte.Rdata")

# FIGURE S2.2
#representation using ggmetanet and diffplot
#ggmetanet (with custom parameters)
ggnet.custom1 = metanetwork::ggnet.default
ggnet.custom1$label = F
ggnet.custom1$edge.alpha = 0.01
ggnet.custom1$max_size = 3
ggnet.custom1$arrow.gap = 0.003
ggnet.custom1$alpha_diff = 0.2
ggnet.custom1$edge.alpha_diff = 0.1

p_metaweb_sp = ggmetanet(g = meta_verte$metaweb,beta = beta,metanetwork = meta_verte,
          ggnet.config = ggnet.custom1,legend = "SBMgroup",
          alpha_per_group = list(groups = diets,
                                 alpha_focal = 0.3,
                                 alpha_hidden = 0.9),
          flip_coords = T)

ggsave(plot = p_metaweb_sp,filename = 'metaweb_sp.tiff',height=2000,width=2500,units = "px")

# at a SBM group level
##################
ggnet.custom2 = ggnet.default
ggnet.custom2$label = T
ggnet.custom2$edge.alpha = 0.5
ggnet.custom2$max_size = 6

ggnet.custom2$label.size = 7
ggnet.custom2$max_size = 17
ggnet.custom2$edge.size = 1
ggnet.custom2$edge.alpha = 0.5
ggnet.custom2$legend.position = "right"
ggnet.custom2$arrow.gap = 0.025

beta = 0.002
# meta_verte = attach_layout(g = meta_verte$metaweb_SBMgroup,
#                   metanetwork = meta_verte,beta = beta)
# save(meta_verte,file = "foodwebs_vs_land_use/data/meta_verte.Rdata")

### FIGURE 2
pGr = ggmetanet(g = meta_verte$metaweb_SBMgroup,
                beta = beta,metanetwork = meta_verte,
          ggnet.config = ggnet.custom2,
          flip_coords = T,
          edge_thrs = 0.8,
          nrep_ly = 1,
          legend = "SBMgroup")

png('metaweb_SBMgroup.png',height=800,width=1200)
print(pGr)
dev.off()

ggnet.custom$label.size=12
ggnet.custom$max_size=17
ggnet.custom$label = T
ggnet.custom$edge.size=1
ggnet.custom$edge.alpha = 0.7
ggnet.custom$edge.alpha_diff=0.7
ggnet.custom$legend.position="right"
ggnet.custom$arrow.gap = 0.025

###FIGURE 4
pDifGr = diff_plot(g1 = meta_verte$inhigh_SBMgroup,
          g2 = meta_verte$inlow_SBMgroup,
          beta = beta,
          metanetwork = meta_verte,
          layout_metaweb = T,
          flip_coords = T,
          edge_thrs = 0.05,
          ggnet.config = ggnet.custom,
          alpha_per_node = list(nodes = dietNames,
                                alpha_focal = 0.4,
                                alpha_hidden = 1))

png('diffplot_marc_group.png',height=800,width=1200)
print(pDifGr)
dev.off()

## FIGURE with Piecharts and ggraph
#(i) building pie chart table
pieChart =
  meta_verte$trophicTable %>% group_by(SBMgroup, Class) %>%
  summarize(nsp = length(unique(species))) %>%
  pivot_wider(names_from = Class, values_from = "nsp") %>%
  ungroup()

pieChart[is.na(pieChart)] = 0
pieChart = cbind(pieChart,x= V(meta_verte$metaweb_SBMgroup)$layout_beta0.002,
                 y = V(meta_verte$metaweb_SBMgroup)$TL,
                 nsp = rowSums(pieChart[2:ncol(pieChart)]))

pieChart$radius = as.vector(scale(log(pieChart$nsp)+1, center=F)*5)
pieChart$radius[c(17,25,29)] = 0.01
options(ggplot2.continuous.colour="viridis")

#(ii) representing the network using scatterpie and ggraph
g = meta_verte$metaweb_SBMgroup
#thresholding on edge
g = delete_edges(g,which(E(g)$weight<0.5))
#tuning parameters
K = 15
R = 0.2

#for scatterpie size legend
seq_pie_size = seq(min(R*pieChart$radius),max(R*pieChart$radius),length = 5)
label_sizes = sapply(seq_pie_size,function(x) min(pieChart$nsp[pieChart$radius>=x/R]))
names(label_sizes) = seq_pie_size

#the plot!
p_pie_graph = ggraph(g, "manual", x = V(g)$layout_beta0.002, y = K*V(g)$TL) +
    geom_edge_link(aes(alpha = weight),arrow = arrow(length = unit(1, 'mm')),
                   end_cap = circle(5, 'mm'),width = 0.15) +
  geom_scatterpie(data = pieChart,
                  aes(x = x,y = K*y,group= SBMgroup, r = R*radius),
                  color = "black", alpha = 0.8,
                  cols = c("Aves", "Mammalia", "Reptilia",
                           "Amphibia", "Basal resources"))+
  scale_fill_manual("Class", values = viridis(5))+
  scale_size_continuous(range = c(0, 0.5)) +
  theme_graph() + 
  geom_scatterpie_legend(seq_pie_size,x = -15, y = 0, labeller = function(x) label_sizes[as.character(x)]) 

p_pie_graph

# FIGURE group layout
## attaching alternative layout 'group-TL-tsne' using SBM level layout
# CAREFUL: you need to attach 'TL-tsne' layout at SBMgroup level before computing group-TL-tsne layout (precomputed here)
group_layout.custom = group_layout.default
group_layout.custom$nbreaks_group = 4
group_layout.custom$group_height = c(1,2,3,4)
group_layout.custom$group_width = c(1,2,3,4)
meta_verte = attach_layout(meta_verte,beta = 0.002,mode = "group-TL-tsne",
                           res = "SBMgroup",group_layout.config = group_layout.custom)

ggnet.custom1 = metanetwork::ggnet.default
ggnet.custom1$label = F
ggnet.custom1$edge.alpha = 0.01
ggnet.custom1$max_size = 3
ggnet.custom1$arrow.gap = 0.003
ggnet.custom1$alpha_diff = 0.2
ggnet.custom1$edge.alpha_diff = 0.1

p_group_metaweb_sp = ggmetanet(g = meta_verte$metaweb,beta = 0.002,metanetwork = meta_verte,
                         mode = "group-TL-tsne",legend = "SBMgroup",
                         ggnet.config = ggnet.custom1,
                         alpha_per_group = list(groups = diets,
                                                alpha_focal = 0.3,
                                                alpha_hidden = 0.9),
                         flip_coords = T)

ggsave(plot = p_group_metaweb_sp,filename = 'group_metaweb_sp.tiff',
       height=2000,width=2500,units = "px")  
  
  
  
### TABLE S2.3
groupis = data.frame(TL=vertex_attr(meta_verte$metaweb_SBMgroup)$TL,
           ab=vertex_attr(meta_verte$metaweb_SBMgroup)$ab,
           gr=vertex_attr(meta_verte$metaweb_SBMgroup)$name)
groupis = groupis[nchar(as.character(groupis$gr))<=2,]
groupis = groupis[order(groupis$TL,decreasing = T),]
#sum(groupis$ab>quantile(groupis$ab,prob=c(.57) ))
#groupis[groupis$ab>quantile(groupis$ab,prob=.57),]
groupis$mostComon = NA
groupis$mostComonClass = NA
groupis$nSpecies = NA
abTot = colSums(meta_verte$abTable)
for(i in 1:dim(groupis)[1]){
  tmp = abTot[rownames(meta_verte$trophicTable)[meta_verte$trophicTable$SBMgroup==groupis$gr[i]]]
  tmp = tmp[order(as.numeric(tmp),decreasing=T)]
  groupis$mostComon[i] = paste(names(tmp)[which.max(tmp)],collapse = ",")
  groupis$nSpecies[i] = length(tmp)
  classes = sapply(names(tmp),function(sp)spTab$Class[spTab$name==sp]) 
  tabi = table(classes)
  groupis$mostComonClass[i] = paste(names(tabi)[tabi==max(tabi)],collapse = ",")
}

groupis=groupis[,!colnames(groupis)%in%c("TL","ab")]
groupis = groupis[c('gr','nSpecies','mostComon','mostComonClass')]
colnames(groupis)[colnames(groupis)%in%c('gr','mostComon','mostComonClass')]=
  c('Group',
    'Most frequent species',
    'Most common class')

save(groupis,file = "groupis")

TMP=data.frame(gr=1:46,onlyBasal=NA,hasBasal=NA)
for(i in 1:46){
  cd = !is.na(spTab$SBMgroup) & spTab$SBMgroup==i
  TMP$onlyBasal[i] = sum(spTab$isBasal[cd])==sum(cd)
  TMP$hasBasal[i] = sum(spTab$isBasal[cd])>0
}
TMP$gr[TMP$onlyBasal]

TMP$gr[TMP$hasBasal]
