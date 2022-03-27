rm(list=ls())
#This script aims at visualising and comparing networks for Botella paper
#using metanetwork package

# setwd("~/Bureau/article_chris_b/")
mainDir ="C:/Users/user/pCloud local/boulot/data/GBIF IUCN tetrapods/resultsRepro" 
setwd(mainDir)

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

# loading data
load("preprocessed_data")

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
spTab_loc = spTab[,c("name","SBMgroup")]
rownames(spTab_loc) = spTab_loc$name
colnames(spTab_loc) = c("species","SBMgroup")
spTab_loc[513,"SBMgroup"] = 16 #assigning NA species to arbitrary group

#adding diets to trophicTable
diets = V(metaweb)$name[ which(!(V(metaweb)$name %in% rownames(spTab_loc)))]
spTab_loc = rbind(spTab_loc,cbind(species = diets,SBMgroup = diets))
rownames(spTab_loc) = spTab_loc$species
#adding diets to abTable (uniform abundance, equals to )
value_diets = min(abTable_loc[which(abTable_loc>0)])
df = data.frame(matrix(rep(value_diets,length(diets)),nrow = 1)) 
df = rbind(df, df[rep(1, nrow(abTable_loc)-1),])
colnames(df) = diets
rownames(df) = rownames(abTable_loc)
abTable_loc = cbind(abTable_loc,df)

# build metanetwork
meta_verte = build_metanet(metaweb = metaweb,abTable = abTable_loc,
                             trophicTable = spTab_loc)
meta_verte = append_agg_nets(meta_verte)
meta_verte = compute_TL(meta_verte)

print_metanet(meta_verte)

#METANETWORK ANALYSIS#
######################

#at a species level
##################

#attaching layout
beta = 0.00005
meta_verte = attach_layout(g = meta_verte$metaweb,
                         metanetwork = meta_verte,beta = beta)

# FIGURE S2.2
#representation using ggmetanet and diffplot
#ggmetanet (with custom parameters)
ggnet.custom = metanetwork::ggnet.default
ggnet.custom$label = F
ggnet.custom$edge.alpha = 0.01
ggnet.custom$max_size = 3
ggnet.custom$arrow.gap = 0.003
ggnet.custom$alpha_diff = 0.2
ggnet.custom$edge.alpha_diff = 0.1
p = ggmetanet(g = meta_verte$metaweb,beta = beta,metanetwork = meta_verte,
          ggnet.config = ggnet.custom,legend = "SBMgroup",
          alpha_per_group = list(groups = diets,
                                 alpha_focal = 0.3,
                                 alpha_hidden = 0.9),
          flip_coords = T)

setwd(mainDir)
png('metaweb_marc_1.png',height=600,width=700)
print(p)
dev.off()

# at a SBM group level
##################
ggnet.custom = ggnet.default
ggnet.custom$label = T
ggnet.custom$edge.alpha = 0.5
ggnet.custom$max_size = 6

ggnet.custom$label.size=7
ggnet.custom$max_size=17
ggnet.custom$edge.size=1
ggnet.custom$edge.alpha = 0.5
ggnet.custom$legend.position="right"
ggnet.custom$arrow.gap = 0.025

beta = 0.002
meta_verte = attach_layout(g = meta_verte$metaweb_SBMgroup,
                  metanetwork = meta_verte,beta = beta)

### FIGURE 2
pGr = ggmetanet(g = meta_verte$metaweb_SBMgroup,
                beta = beta,metanetwork = meta_verte,
          ggnet.config = ggnet.custom,
          flip_coords = T,
          edge_thrs = 0.8,
          nrep_ly = 1)

png('metaweb_marc_group.png',height=800,width=1200)
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

### TABLE S2.3
groupis=data.frame(TL=vertex_attr(meta_verte$metaweb_SBMgroup)$TL,
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
