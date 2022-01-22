# ADAPT directories 
repoDir = "C:/Users/user/pCloud local/boulot/Github/tetrapods_trophicNet_VS_anthropization/"
setwd(repoDir)
source('functions_to_source.R')

occDir = "C:/Users/user/pCloud local/boulot/data/GBIF IUCN tetrapods/"
saveDir = paste(occDir,'resultsRepro/',sep="")

setwd(saveDir)
load(file="preprocessed_data")
load(file = 'TrophicNetworksList')

#####
# species trophic statistics from metaweb
#####

### Get metaweb position statistics
diag(metaweb_AdjAdult)=0 
Adj_tmp = metaweb_Adj;diag(Adj_tmp)=0
# Trophic level MacKay
G = graph_from_adjacency_matrix(t(Adj_tmp),mode = "directed")
tls = trophicLevelMackay(G)
stats = data.frame(name=rownames(metaweb_Adj),tl=tls) 

# We consider basal and prey among tetrapods only  
tetra = !colnames(metaweb_Adj)%in%dietNames
Adj_TetraOnly = metaweb_Adj[tetra,tetra]
stats = stats[!stats$name%in%dietNames,]

stats$isBasal = sapply(1:dim(Adj_TetraOnly)[1],function(i)sum(Adj_TetraOnly[i,])==0)
# Count number of preys (all life stages)
stats$nPrey = sapply(1:dim(Adj_TetraOnly)[1],function(i)sum(Adj_TetraOnly[i,]))
# Count number of predator on the adult only
stats$nPred = sapply(1:dim(metaweb_AdjAdult)[1],function(i)sum(metaweb_AdjAdult[,i]))
stats$isTop = stats$tl>2.262

tlBreaks = c(0,max(stats$tl[stats$isBasal]),
             min(stats$tl[stats$isTop])-.001,max(stats$tl))
stats$tlforOmni= base::cut(stats$tl,
                           breaks=tlBreaks)

#forLine = expand.grid(y=c(0,100),x=tlBreaks)
#ggplot()+geom_histogram(data=stats,aes(x=tl),bins=70)+
#  geom_line(data=forLine,aes(x=x,y=y,group=x))

# Omnivory levels
stats$omniLev = NA;stats$isOmni=NA;
for(i in 1:dim(stats)[1]){
  cd = metaweb_Adj[stats$name[i],]>0
  if(sum(cd)>1){
    stats$omniLev[i] = sd(tls[cd])
  }else{stats$omniLev[i] = 0 # For basals and those with only one prey
  }
  prey_tmp = colnames(metaweb_Adj)[cd]
  if(length(prey_tmp)>2){
    # We look consider omnivore species feeding on all regrouped trophic levels 
    stats$isOmni[i] = length(unique(stats$tlforOmni[stats$name%in%prey_tmp]))>1
  }else{stats$isOmni[i] = F}
}


#rownames(metaweb_Adj)[stats$isTop]
spTab = spTab[,!colnames(spTab)%in%c("tl","omniLev","isOmni","isBasal","nPrey","nPred","isTop",'tlforOmni')]

spTab = merge(spTab,stats,by="name",all.x=T)

#####
# Compute metrics
#####

metrics = unlist(ebdCompo[1:6])
ebd = data.frame(cell=as.numeric(names(gList)))
for(met in as.character(metrics)){eval(parse(text=paste('ebd$',met,'=NA',sep="")))}

time.in.mins= function(time){ 60*as.numeric(format(time, format='%H'))+
    as.numeric(format(time, format='%M'))+
    as.numeric(format(time, format='%S'))/60}

### Top/Bascal proportions, basal to top ShortestPaths statistics
deb = time.in.mins(Sys.time())
for(i in 1:length(gList)){
  g = gList[[i]]
  
  assemblage = igraph::get.vertex.attribute(g,name = "name")
  nSp = length(assemblage)
  adj = as.matrix(as_adj(g))
  adj_adultWeb = t(metaweb_AdjAdult[assemblage,assemblage])
  
  # proportion of apex (metaweb)
  ebd$pApexMeta[i] = sum(spTab$isTop[spTab$name%in%assemblage])/ nSp
  # proportion of basal (metaweb)
  ebd$pBasalMeta[i] = sum(spTab$isBasal[spTab$name%in%assemblage])/ nSp
  # proportion of basal (local)
  basal = colSums(adj)==0
  ebd$pBasal[i] = sum(basal)/ nSp
  # dirCON!!!!
  ebd$dirCon[i] = sum(adj)/(nSp^2)
  
  # Omnivory metrics (corrected for basals)
  cd = spTab$name%in%assemblage & !spTab$isBasal & !spTab$isTop
  ebd$omniLev[i] = mean(spTab$omniLev[cd],na.rm=T)
  ebd$omniProp[i] = sum(spTab$isOmni[cd],na.rm = T)/sum(cd)
  
  apex = rowSums(adj_adultWeb)==0 & colSums(adj)>0
  # Matrix of shortest-paths lengths between basal and top species 
  D = igraph::distances(g,mode="out")
  if(sum(apex)>0 & sum(basal)>0){
    basalToApex = as.vector(D[as.logical(basal),as.logical(apex),drop=F] )
    pathsLengths = basalToApex[is.finite(basalToApex)]
    ebd[i,c('maxPath','meanPath','sdPath')] = c(max(pathsLengths),mean(pathsLengths),sd(pathsLengths))
  }else{ebd[i,c('maxPath','meanPath','sdPath')]=rep(0,3)}

  # Modularity 
  adj_undir = 1. * matrix(as.vector(adj + t(adj))>0,dim(adj)[1],dim(adj)[1])
  g_undir = graph_from_adjacency_matrix(adj_undir, mode = "undirected")
  if(T){
    degs = igraph::degree(g_undir)
    membership = igraph::components(g_undir,mode = "weak")$membership
    globMembership = as.character(membership)
    for(m in unique(membership)){
      members = which(membership==m)
      if(length(members)>1){
        g_sub = igraph::induced.subgraph(g_undir,vids=members)
        spinGlassClust = cluster_spinglass(g_sub)
        globMembership[members] = paste(m,spinGlassClust$membership)
      }else{globMembership[members] = as.character(m)}
    }
    nLinks = sum(adj_undir[lower.tri(adj_undir)])
    comMatrix = 0. * adj_undir
    for(m in unique(globMembership)){
      comMatrix[which(globMembership==m),which(globMembership==m)] = 1
    }
    expectedRandomLinkProbaMatrix = matrix(degs,length(degs),1) %*% matrix(degs,1,length(degs)) / (2*nLinks-1)
    termsMat = comMatrix * (adj_undir - expectedRandomLinkProbaMatrix)
  }
  ebd$modul[i] = sum(termsMat[lower.tri(termsMat)])/(2*nLinks)
  # meanShortDist 
  D = igraph::distances(g_undir)
  dists = as.vector(D)
  dists = dists[is.finite(dists) & dists!=0]
  ebd$meanShortDist[i] = mean(dists)
  
  if(i/30==round(i/30)){
    setwd(saveDir)
    write.table(ebd,'fullEbd.csv',sep=";",row.names=F,col.names=T)
    dur = time.in.mins(Sys.time())-deb
    flush.console()
    cat('\r Processed...',round(1000*i/length(gList))/10,'% in ',dur,' minutes')
  }
}

setwd(saveDir)
write.table(ebd,'fullEbd.csv',sep=";",row.names=F,col.names=T)


#####
# Compute fragmentation indices
#####

setwd(saveDir)
ebd=read.csv('fullEbd.csv',sep=";",header=T,stringsAsFactors = F)

MAT=as.matrix(rCells)
R = as.matrix(rLandUse)

rayon=4
ebd$patchAntiArea=NA
ebd$proxToBorder=NA
ebd$divLandUse=NA

D = matrix(NA,2*rayon+1,2*rayon+1)
for(k in 1:dim(D)[1]){
  D[k,] = sqrt( (k-(rayon+1))^2+(1:dim(D)[2]-(rayon+1))^2 )
}

for(i in 1:dim(ebd)[1]){
  coos=c(NA,NA);
  pos = which(MAT==as.numeric(ebd$cell[i]))
  coos[1] = pos%%dim(MAT)[1]
  coos[2] = 1+(pos%/%dim(MAT)[1])
  mat = R[(coos[1]-rayon):(coos[1]+rayon),(coos[2]-rayon):(coos[2]+rayon)]
  
  # Area of land use (land system + intensity) patch (of the land use of central pixel)
  extr = mat;extr[] = as.numeric(mat[]==mat[rayon+1,rayon+1])
  extr = ConnCompLabel(extr)
  targetCompo = extr[rayon+1,rayon+1]
  ebd$patchAntiArea[i] = -sum(as.vector(extr)==targetCompo,na.rm = T)
  # Distance to edge 
  if(T){
    extr[] = as.numeric(!is.na(as.vector(extr)) & as.vector(extr)==targetCompo)
    
    m1 = convolution(extr,
                     matrix(c(0,0,0,
                              0,1,-1,
                              0,0,0),3,3))
    m2 = convolution(extr,
                     matrix(c(0,0,0,
                              -1,1,0,
                              0,0,0),3,3))
    m3 = convolution(extr,
                     matrix(c(0,-1,0,
                              0,1,0,
                              0,0,0),3,3))
    m4 = convolution(extr,
                     matrix(c(0,0,0,
                              0,1,0,
                              0,-1,0),3,3))
    border = extr
    border[] = as.numeric(as.vector(m1+m2+m3+m4)>0)
    bordVec = as.vector(border)==1
    if(sum(bordVec)>0){
      distsToBorder = D[bordVec]  
      ebd$proxToBorder[i] = -distsToBorder[which.min(distsToBorder)][1]
    }else{
      ebd$proxToBorder[i] = -(rayon+1)
    }
  }
  
  # Diversity of surrounding land uses
  surroundingLandUse = unique(as.vector(mat[(rayon:(rayon+2)),(rayon:(rayon+2))]))
  surroundingLandUse = surroundingLandUse[!is.na(surroundingLandUse)]
  ebd$divLandUse[i] = length(surroundingLandUse)
  #print(ebd[i,c('patchAntiArea','proxToBorder','divLandUse')])
  
  if(i/50==round(i/50)){
    flush.console()
    cat('\r Process...',round(i*1000/dim(ebd)[1])/10,'%')
  }
}

setwd(saveDir)
write.table(ebd,'fullEbd_frag.csv',sep=";",row.names=F,col.names=T)

#####
# Metrics correlations
#####

setwd(saveDir)
met = read.csv('fullEbd.csv',sep=";",header=T,stringsAsFactors = F)

rownames(met) = met$cell

met = met[complete.cases(met),]
ordered = unlist(ebdCompo[1:6])
formu = as.formula(paste("~",paste(ordered,collapse="+")))
samp = sample(1:dim(met)[1],650)
setwd(saveDir)
png("metrics_pairs.png",height=1600,width=1600)
pairs( formu ,data=met[samp,] , 
       diag.panel = panel.hist, 
       upper.panel = NULL , 
       labels= ordered,
       pch=3 ,
       col=  rep('black',length(samp)) ,
       cex.labels = 1.6 , 
       cex.axis = 2)
dev.off()

# Correlations across metrics
mcor = matrix(NA,length(ordered),length(ordered))
for(i in 1:length(ordered)){
  for(j in 1:length(ordered)){
    mcor[i,j] = cor(met[,ordered[i]] , met[,ordered[j]] )
  }
}
rownames(mcor) = ordered
colnames(mcor) = ordered

setwd(saveDir)
png('metrics_correlations.png',height=1000,width=1200)
print(corrplot(mcor, type="upper", tl.col="black", tl.srt=45, diag=F,cl.cex=1.5,cl.lim=c(-1,1),addCoef.col = "grey80",number.cex=1))
dev.off()

#####
# Design matrix
#####

setwd(saveDir)
ebd = read.csv('fullEbd_frag.csv',sep=";",header=T,stringsAsFactors = F)
rownames(ebd) = ebd$cell
minSamp = 29

# Filter rows with NA in embedding
Y = ebd[,colnames(ebd)!="cell"]
Y = Y[as.character(cells$cell),]
hasY = complete.cases(Y[,unlist(ebdCompo)])
Y = Y[hasY,unlist(ebdCompo)]
design = cells[hasY,]

# PLOT number of cells per region and land use
if(F){
  tmpRegions = regions
  tmpRegions$xCoo = c(4,NA,NA,1,NA,5,2,NA,3,NA,NA,NA)
  keptLandUses$yCoo = 17:1
  toPlot = as.data.frame(expand.grid(xCoo=1:9,yCoo=1:17))
  toPlot = merge(toPlot,tmpRegions[,c('xCoo','region')],by="xCoo")
  toPlot = merge(toPlot,keptLandUses[,c('yCoo','code')],by="yCoo")
  toPlot$nCells = 0
  
  cases = aggregate(list(count=rep(1,dim(design)[1])),
                    by=list(landGroup=design$landGroup,landUse=design$landUse,
                            region=design$region,
                            intens=design$intens),sum)
  for(i in 1:dim(cases)[1]){
    toPlot$nCells[toPlot$region==cases$region[i] & toPlot$code==cases$landUse[i]] = cases$count[i]
  }
  
  toPlot$text = as.character(toPlot$nCells)
  toPlot$colorGroup = sapply(1:dim(toPlot)[1],function(i){
    if(toPlot$nCells[i]==0){y="red"
    }else if(toPlot$nCells[i]<10){y="orange"
    }else if(toPlot$nCells[i]<30){y="yellow"
    }else{y="olivedrab2"}
    return(y)})
  toPlot = toPlot[,c('yCoo','xCoo','text','colorGroup')] 
  
  ColTitles = data.frame(yCoo=18.5,
                         xCoo=1:sum(!is.na(tmpRegions$xCoo)),
                         text=sapply(1:sum(!is.na(tmpRegions$xCoo)),function(x){
                           tmpRegions$name[!is.na(tmpRegions$xCoo) & tmpRegions$xCoo==x]}),
    colorGroup="white")
  toPlot = rbind(toPlot,ColTitles)
  RowTitles = data.frame(yCoo=1:17,xCoo=-.5,text="",
                         colorGroup= sapply(1:17,function(x){keptLandUses$color[keptLandUses$yCoo==x]}) )
  toPlot = rbind(toPlot,RowTitles)
  colorLevels = c('red','orange','yellow','olivedrab2',as.character(keptLandUses$color),'white')
  toPlot$colorGroup = factor(toPlot$colorGroup,levels=colorLevels)
  p = ggplot(toPlot,aes(x=xCoo,y=yCoo,fill=colorGroup,label=text))+geom_tile()+geom_text(size=7)
  p = p + scale_fill_manual(values=colorLevels)
  p = p +theme_bw()+xlab('Bioclimatic region')+ylab('Land system and use intensity')
  p = p + theme(legend.position = "none",axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank(),axis.title = element_text(size=25))
  
  setwd(saveDir)
  png('nCellsPerLandGroup_taxoCoverSup70_richSup20.png',height=900,width=1300)
  print(p)
  dev.off()
}

# Filter rows belonging to scarcely sampled cases
toKeep = rep(T,dim(design)[1])
for(reg in unique(design$region)){
  for(lsy in unique(design$landSys)){
    cd1 = design$region==reg & design$landSys==lsy 
    underLow = sum(cd1 & design$intens=="low")
    if(underLow>minSamp){
      for(inte in c('mid','high')){
        nbCells = sum(cd1 & design$intens==inte)
        if(nbCells<minSamp){toKeep[cd1 & design$intens==inte]=F}
      }
    }else{toKeep[cd1] = F}
  }
}

design = design[toKeep,]
Y = Y[toKeep,]

design$intens = factor(design$intens , levels=c('low','mid','high'))
design$intens = C(design$intens,contr.treatment)
design$landSys = factor(design$landSys, levels=c("forest",as.character(unique(design$landSys)[unique(design$landSys)!="forest"])))
design$landSys = C(design$landSys, contr.treatment)
tmpRegs = as.numeric(as.character(design$region))
design$region = factor(tmpRegs,
    levels=c(7,unique(tmpRegs)[unique(tmpRegs)!=7]))
design$region = C(design$region , contr.treatment)
preX = design[,c("region","landSys",'intens')]
formu=as.formula(paste('~ ',
  paste(c('region','landSys','intens'),
  collapse = "/")))
X = model.matrix(object=formu,
                 data = preX)
X = X[,colSums(X)>0]

setwd(saveDir)
save(X,Y,file = 'MultiRegMatrices')
save(design,file = 'study_design_tab')

#####
# Multivariate centroid estimates vs  bioclimatic region, land use and intensity
#####

setwd(saveDir)
load(file = 'MultiRegMatrices')
load(file = 'study_design_tab')

### Compute centroids along metrics per land group
XX= t(X)%*%X
XXinv = solve(XX)
B = XXinv %*% t(X) %*% as.matrix(Y)
rownames(B) = colnames(X)
#Y_ = X %*% B 

# Change effects encoding name
rownames(B)[rownames(B)=="(Intercept)"] = "(Intercept)_region7:landSysforest:intens1"
contr = attr(design$region,which = "contrasts")
cd = regexpr('region',rownames(B))>0 & nchar(rownames(B))==7
for(elem in which(cd)){
  reg = rownames(contr)[contr[,substr(rownames(B)[elem],7,7)]==1]
  rownames(B)[elem] = paste(substr(rownames(B)[elem],1,6),reg,sep="")
}
contr = attr(design$landSys,which = "contrasts")
cd = regexpr(':landSys',rownames(B))>0 & nchar(rownames(B))==16
for(elem in which(cd)){
  ls = rownames(contr)[contr[,substr(rownames(B)[elem],16,16)]==1]
  rownames(B)[elem] = paste(substr(rownames(B)[elem],1,15),ls,sep="")
}
contr = attr(design$intens,which = "contrasts")
cd = regexpr(':intens',rownames(B))>0 & regexpr('Intercept',rownames(B))<=0
for(elem in which(cd)){
  lenCh = nchar(rownames(B)[elem])
  ls = rownames(contr)[contr[,substr(rownames(B)[elem],lenCh,lenCh)]==1]
  rownames(B)[elem] = paste(substr(rownames(B)[elem],1,lenCh-1),ls,sep="")
}

setwd(saveDir)
saveRDS(B,file="B_pointEstim")

#####
# Figure 3: Relative deviation related to higher intensity per metric 
#####

setwd(saveDir)
B=readRDS(file="B_pointEstim")
load(file = 'MultiRegMatrices')
load(file = 'study_design_tab')
cd = design$intens!="low"
data = design[cd,]
Y = Y[cd,]

res = as.data.frame(expand.grid(metric=unlist(ebdCompo[1:6]),landGroup=unique(data$landGroup),dev=NA))
res$landGroup = as.character(res$landGroup)
res = merge(res,unique(data[,c('landGroup','intens','landUse','landSys','region')]),by="landGroup",all.x=T)
for(i in 1:dim(res)[1]){
  for(j in 1:length(ebdCompo)){
    if(res$metric[i]%in%ebdCompo[[j]]){
      res$ebd[i] = names(ebdCompo)[j]
    }
  } 
}
for(i in 1:dim(res)[1]){
  encod = paste('region',res$region[i],':landSys',res$landSys[i],':intens',res$intens[i],sep="")
  if(encod %in% rownames(B)){
    RelDev =  B[encod,res$metric[i]]
    res$dev[i] = RelDev / IQR(Y[,res$metric[i]],na.rm = T)
  }
}
res$agree='mean region and land system'
res$hypo=NA
res$metric = paste('\n ',res$metric,sep="")
res$metric=factor(res$metric,levels=paste('\n ',as.character(unlist(ebdCompo[1:6])),sep=""))

resSum = aggregate(list(dev=res$dev),by=list(metric=res$metric,ebd=res$ebd,intens=res$intens),mean)
resSum$hypo = sapply(1:dim(resSum)[1],function(i){
  if(resSum$ebd[i]=="basal"){"-"}
  else if(resSum$ebd[i]=="apex"){"-"}
  else if(resSum$ebd[i]=="connectance"){"+"}
  else if(resSum$ebd[i]=="omnivory"){"+"}
  else if(resSum$ebd[i]=="vertiChains"){"-"}
  else if(resSum$ebd[i]=="compart"){"-"}
  else if(resSum$ebd[i]=="fragment"){"+"}
})
resSum$agree = as.character((resSum$hypo=="+" & resSum$dev>0) | (resSum$hypo=="-" & resSum$dev<0))

if(F){
  toPlot = rbind(resSum[,c('metric','ebd','intens','dev','hypo','agree')],
                 res[c('metric','ebd','intens','dev','hypo','agree')])
  toPlot$agree = factor(toPlot$agree,levels=c("TRUE","FALSE","mean region and land system"))
  
  p = ggplot(toPlot[toPlot$intens=="high",],aes(x=metric,y=dev,colour=agree,size=agree))+geom_point()
  p = p+scale_colour_manual(values=c('green','red','grey'))+scale_size_manual(values=c(3,3,1))
  p = p +xlab('Network metric')+ylab('Relative deviation (from low land use intensity)')+theme_bw()
  print(p)
}

resSum$agree = factor(resSum$agree,levels=c('TRUE','FALSE'))

plotList = lapply(c('mid','high'),function(intensity){
  p = ggplot()+geom_violin(data=res[res$intens==intensity,],aes(x=metric,y=dev),colour="grey")+
    geom_point(data=res[res$intens==intensity,],aes(x=metric,y=dev),colour='grey')+
    geom_col(data=resSum[resSum$intens==intensity,],aes(x=metric,y=dev,fill=agree),alpha=.4)+
    scale_fill_manual(name="Agreement with hypothesis",values=c('green','red'))+
    ylab(paste('Relative deviation \n ',intensity,' minus low land use intensity',sep=""))+
    theme_bw()+theme(axis.title.x = element_blank(),text =element_text(size=20),axis.text.x = element_text(size=20,angle=20))
  return(p)
})

setwd(saveDir)
png('deviation_related_to_intensity.png',height=750,width=1000)
multiplot(plots=list(plotList[[2]]),cols = 1)
dev.off()

#####
# Compute pValues
#####

# real data
setwd(saveDir)
load(file = 'MultiRegMatrices')
load(file = 'study_design_tab')
data= design
Y = as.matrix(Y)

if(F){
  ### Are there significant differences across land systems AND regions? 
  # compair pairs only
  formula = Y~regAndLandSys;test = c(0, 0, 0, 1)
  alpha = 0.05;factors.and.variables = T
  formula = Formula::Formula(formula)
  frame = model.frame(formula, data = data)
  groupvar.location = length(frame[1, ])
  groupvar = names(frame)[groupvar.location]
  vars = names(frame)[1:(groupvar.location - 1)]
  if (!is.factor(frame[, groupvar])) {
    frame[, groupvar] <- factor(frame[, groupvar])
  }
  levels = levels(frame[, groupvar])
  o = order(frame[, groupvar])
  frame = frame[o, ]
  p = length(vars)
  a = length(levels(frame[, groupvar]))
  k = a - 1
  n=1
  step2subsets = vector("list", a*(a-1)/2)
  for(i in 1:(a-1)){
    for(j in (i+1):a){
      subsetframe <- subset(frame, frame[, groupvar] %in% c(levels[c(i,j)]))
      subsetframe <- droplevels(subsetframe)
      groupvarsub = names(subsetframe)[groupvar.location]
      base <- nonpartest(formula = formula,
                         data= subsetframe, 
                         tests = test,
                         plots = F,permtest = F)
      if (test[1] == 1) {
        testpval = base$results$`P-value`[1]
      }
      if (test[2] == 1) {
        testpval = base$results$`P-value`[2]
      }
      if (test[3] == 1) {
        testpval = base$results$`P-value`[3]
      }
      if (test[4] == 1) {
        testpval = base$results$`P-value`[4]
      }
      if (testpval * a/k < alpha) {
        step2subsets[[n]] = levels(subsetframe[, 
                                               groupvarsub])
      }
      else{step2subsets[[n]] = NA}
      if(n/10==round(n/10)){
        flush.console()
        cat('\r Processed...',round(1000*n*2/(a*(a-1)))/10,'%')
      }
      n=n+1
    }
  }
  step2subsets = step2subsets[!is.na(step2subsets)]
  
  setwd(saveDir)
  saveRDS(step2subsets,'regionsAndLandSys_pairs_tests')
  
  # Are there significant differences across land systems? -> YES across all pairs
  # Are there significant differences across regions? -> YES across all pairs
  # (without accounting for other variabilities)
  for(factorToTest in c("region",'landSys')){
    print(paste('differences across ',factorToTest))
    dataTmp = data
    dataTmp[,factorToTest] = factor(dataTmp[,factorToTest])
    Y = ebd[ebd$cell%in%dataTmp$cell,!colnames(ebd)=="cell"]
    Y = as.matrix(Y[as.character(dataTmp$cell),])
    cd = complete.cases(Y);Y = Y[cd,];dataTmp = dataTmp[cd,];
    tmp = ssnonpartest( as.formula(paste('Y~',factorToTest)), 
                        data = dataTmp,  
                        test = c(0, 0, 0, 1), alpha = 0.05, 
                        factors.and.variables = T)
  }
  
}

# Are there differences across intensity? Inside each region AND landSys
results = unique(data[,c('region','landSys','intens')])
results=results[results$intens!="low",]
results$testUnit = paste(results$region,'_',results$landSys,'_',results$intens,sep="")
expEbds = as.data.frame(expand.grid(testUnit=results$testUnit,ebd=names(ebdCompo)))
results = merge(expEbds,results,by="testUnit",all.x=T)
results = results[,colnames(results)!="testUnit"]
# We control the probability of 
# at least one fake rejection among all embeddings
alpha_ = 1-.95^(1/length(ebdCompo))
nGroups = dim(unique(results[,c('region','landSys','ebd')]))[1]  
results$sigDiffFromLow = NA
n=1
for(reg in unique(results$region)){
  for(ld in unique(results$landSys)){
    cdRes1 = results$region==reg & results$landSys==ld
    cd = data$landSys==ld & data$region==reg
    if(sum(cdRes1)>0 & sum(cd)>1 & length(unique(data$intens[cd]))>1){
      X = data[cd,]
      for(eb in names(ebdCompo)){
        Ytmp = Y[cd,ebdCompo[eb][[1]]]
        formula = Ytmp~intens
        res = ssnonpartest2(formula,
                            data = X,
                            test = c(0, 0, 0, 1),
                            alpha = alpha_, 
                            factors.and.variables = T)
        for(contrasto in results$intens[cdRes1 & results$ebd==eb]){
          extr = sapply(res,function(el) length(el)==2 & "low"%in%el & contrasto%in%el)
          if(sum(extr==1)){isSig = T}else{isSig = F}
          results$sigDiffFromLow[cdRes1 & results$ebd==eb & results$intens==contrasto] = isSig 
        }
        flush.console()
        cat('\r Processed ',round(1000*n/nGroups)/10,'%')
        n = n+1
      }
      setwd(saveDir)
      write.table(results,'test_intens_per_RegAndLandSys.csv',sep=";",row.names=F,col.names=T)
    }
  }
}

#####
# Figures S4.4 to 4 & S6 (frag) App: Plot tables of effects 
#####

setwd(saveDir)
res = read.csv('test_intens_per_RegAndLandSys.csv',sep=";",header=T,stringsAsFactors = F)
B=readRDS(file="B_pointEstim")
load(file = 'MultiRegMatrices')
load(file = 'study_design_tab')
# remove BLACK and STEPPIC as they lack 
# levels for our comparisons  
cd = !design$region%in%c(11,12,5) & design$intens!="low"
data = design[cd,]
Y = Y[cd,]

### data.frame pour geom_tile
tmpLU = keptLandUses[keptLandUses$intens!="low",]
res$y = sapply(1:dim(res)[1],function(i){dim(tmpLU)[1]-which(tmpLU$intens==res$intens[i] & tmpLU$landSys==res$landSys[i])})
ordReg = c(4,7,9,1,6)
res$x = sapply(1:dim(res)[1],
               function(i){which(ordReg==res$region[i])})
colnames(res)[colnames(res)=="sigDiffFromLow"] = "coloriage"
res$text = NA
res$textColor = "black"
#res$textSize = 4.3
res$textSize = 4.8
res$hypo = sapply(1:dim(res)[1],function(i){
  if(res$ebd[i]=="basal"){"-"}
  else if(res$ebd[i]=="apex"){"-"}
  else if(res$ebd[i]=="connectance"){"+"}
  else if(res$ebd[i]=="omnivory"){"+"}
  else if(res$ebd[i]=="vertiChains"){"-"}
  else if(res$ebd[i]=="compart"){"-"}
  else if(res$ebd[i]=="fragment"){"+"}
})
# row names
initTab = tmpLU[,c('short','color')]
colnames(initTab)=c('text',"coloriage")
initTab$y = (dim(initTab)[1]-1):0
initTab$x = 0
initTab$text=as.character(initTab$text);initTab$coloriage=as.character(initTab$coloriage)
initTab$textSize= 4.2
initTab$textColor = "black"
# col names
initTab2 = data.frame(text=sapply(ordReg,function(reg)regions$name[regions$region==reg]),
                      coloriage=NA,
                      y=max(res$y,na.rm=T)+1,
                      x=1:length(ordReg))
initTab2$text=as.character(initTab2$text)
initTab2$textSize= 6.5
initTab2$textColor = "black"

for(ebd in names(ebdCompo)){
  tmp = res[res$ebd==ebd,c('text','coloriage','y','x','textSize','textColor','hypo','region','landSys','intens')]
  # encode effect of intens. 
  # on all embedding dimensions
  for(i in 1:dim(tmp)[1]){
    encod = paste('region',tmp$region[i],':landSys',tmp$landSys[i],':intens',tmp$intens[i],sep="")
    if(encod %in% rownames(B)){
      coefsVar =  B[encod,ebdCompo[ebd][[1]]]
      names(coefsVar) = ebdCompo[ebd][[1]]
      if(ebd!="fragment"){
        for(k in 1:length(coefsVar)){
          coefsVar[k] = coefsVar[k]/IQR(Y[,names(coefsVar)[k]],na.rm = T)          
        }
      }
      vals = as.character(round(coefsVar*1000)/1000)
      vals = sapply(vals,function(el){
        if(regexpr('-',el)<=0){paste('+',el,sep="")
        }else{el}})
      if(length(vals)>=3){vals[3]= paste(vals[3],'\n',sep="")}
      tmp$text[i] = paste(vals,collapse=";")
    }
    else{print(paste('Missing effect for row ',i))}
  }
  
  # encode "fit expectation" or "contradict exp." through text color
  for(i in 1:dim(tmp)[1]){
    nPlus = sum( gregexpr("\\+", tmp$text[i])[[1]]!=(-1) )
    nMoins = sum( gregexpr("\\-", tmp$text[i])[[1]]!=(-1) )
    nCoos = length(ebdCompo[names(ebdCompo)==ebd][[1]])
    if((nPlus==nCoos & tmp$hypo[i]=="+" & tmp$coloriage[i]==T) | (nMoins==nCoos & tmp$hypo[i]=="-" & tmp$coloriage[i]==T)){
      # Fit hypothesis
      #tmp$textSize[i] = 5.7
      tmp$textColor[i] = "green4"
    }else if((nPlus==nCoos & tmp$hypo[i]=="-" & tmp$coloriage[i]==T) | (nMoins==nCoos & tmp$hypo[i]=="+" & tmp$coloriage[i]==T)){
      # Meet hypothesis
      tmp$textColor[i] = "red1"
    } 
  }
  
  # 
  toPlot = rbind(initTab,tmp[,1:6],initTab2)
  toPlot$coloriage = factor(toPlot$coloriage,levels=c('FALSE','TRUE',initTab$coloriage))
  p = ggplot(toPlot,aes(x=x,y=y,fill=coloriage,label=text))
  p = p+geom_tile(color="black")+geom_text(aes(size=factor(toPlot$textSize),colour=factor(toPlot$textColor)))
  p = p + scale_fill_manual(values=c('grey','deepskyblue1',as.character(tmpLU$color) ))
  p = p + scale_size_manual(values= as.numeric(levels(factor(toPlot$textSize))) )
  p = p + scale_colour_manual(values=c(levels(factor(toPlot$textColor))))
  p = p+theme(legend.position="none",
              panel.background =element_rect(fill="white"),
              panel.grid = element_line(color="white"),
              axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank())
  p = p + ylab('Land system and use intensity') + xlab('bioclimatic region')
  
  setwd(saveDir)
  png(paste('effect_intensity_on_',ebd,'.png',sep=""),height=600,width=1100)
  print(p)
  dev.off()
  
  res[res$ebd==ebd,c('text','coloriage','y','x','textSize','textColor','hypo','region','landSys','intens')] = tmp
}

setwd(saveDir)
write.table(res,'deviations_vs_intensity_Table.csv',sep=";",row.names=F,col.names=T)

######
# Table 3
######

setwd(saveDir)
res=read.csv('deviations_vs_intensity_Table.csv',sep=";",header=T,stringsAsFactors = F)

tab=data.frame(ebd=names(ebdCompo)[1:6],nSig=NA,hypoT=NA,hypoF=NA)
for(i in 1:dim(tab)[1]){
  tab$nSig[i] = sum(res$coloriage[res$ebd==tab$ebd[i]])
  tab$hypoT[i] = round(1000*sum(res$textColor[res$ebd==tab$ebd[i]]=="green4")/ tab$nSig[i])/10
  tab$hypoF[i] = round(1000*sum(res$textColor[res$ebd==tab$ebd[i]]=="red1")/ tab$nSig[i])/10
}
print(tab)

#####
# Figure 13 App: SP length distributions high vs low intensity
#####

setwd(saveDir)
load(file="preprocessed_data")
load(file = 'TrophicNetworksList')

### Top/Bascal proportions, basal to top ShortestPaths statistics
maxLen = 30
lenVals = c(Inf,0:maxLen)
for(i in 1:length(gList)){
  tmp = data.frame(cell=as.numeric(names(gList)[[i]]))
  
  assemblage = igraph::get.vertex.attribute(gList[[i]],name = "name")
  nSp = length(assemblage)
  adj = as.matrix(as_adj(gList[[i]]))
  adj_adultWeb = t(metaweb_AdjAdult[assemblage,assemblage])
  
  # local basals 
  basal = colSums(adj)==0
  # local apex
  apex = rowSums(adj_adultWeb)==0 & colSums(adj)>0
  
  # Matrix of shortest-paths lengths between basal and top species 
  D = igraph::distances(gList[[i]],mode="out")
  if(sum(apex)>0 & sum(basal)>0){
    basalToApex = as.vector(D[as.logical(basal),as.logical(apex),drop=F] )
    counts = table(c(basalToApex,lenVals))-rep(1,length(lenVals))
    counts = counts[as.character(lenVals)]
    #pathsLengths = basalToApex[is.finite(basalToApex)]
    tmp[,paste('len_',lenVals,sep="")] = counts
  }else{tmp[,paste('len_',lenVals,sep="")]=rep(0,length(lenVals))}
  
  if(i==1){ebd = tmp
  }else{ebd = rbind(ebd,tmp)}
  
  if(i/30==round(i/30)){
    setwd(saveDir)
    write.table(ebd,'test_SP_ebd.csv',sep=";",row.names=F,col.names=T)
    flush.console()
    cat('\r Process...',round(1000*i/length(gList))/10,'%')
  }
}
setwd(saveDir)
write.table(ebd,'test_SP_ebd.csv',sep=";",row.names=F,col.names=T)

### PLOT histogram of SP lengths in low vs high intensity regions
setwd(saveDir)
ebd=read.csv('test_SP_ebd.csv',sep=";",header=T)
load(file = 'MultiRegMatrices')
load(file = 'study_design_tab')

rownames(ebd)= ebd$cell
#ebd = ebd[rownames(Y),!colnames(ebd)%in%c("cell","len_Inf")]
ebd = ebd[rownames(Y),!colnames(ebd)%in%c("cell")]

groupFrom = 5
eval(parse(
  text=paste('ebd$len_',groupFrom,'andMore = sapply(1:dim(ebd)[1],function(i)sum(ebd[i,paste("len_",',groupFrom,':30,sep="")],na.rm=T) )',sep="")))

ebd = ebd[,!colnames(ebd)%in%paste('len_',groupFrom:30,sep="")]
ebd = ebd/rowSums(ebd)
colsToKeep = colnames(ebd)[colSums(ebd,na.rm=T)>0] 
toPlot = data.frame(cell=rownames(ebd),
                    count=ebd[,colsToKeep[1]],
                    len=colsToKeep[1])
for(col in colsToKeep){
  toPlot = rbind(toPlot,
                 data.frame(cell=rownames(ebd),
                            count=ebd[,col],
                            len=col))
}

toPlot = merge(toPlot,design[,c('cell','landGroup','intens','region','landSys')],by="cell",all.x=T)
grpToRemove = c('4_31','4_32','7_31','9_31',
                '7_51','9_51','6_51','1_51',
                '1_61',
                '7_76','9_76','1_76','6_76')
tmp = toPlot[toPlot$intens%in%c('low','high') & !toPlot$landGroup%in%grpToRemove,] 
# Mean per landGroup
tmp = aggregate(list(prop=tmp$count),
                by=list(len=tmp$len,
                        intens=tmp$intens,
                        landGroup=tmp$landGroup,
                        region=tmp$region,
                        landSys=tmp$landSys),
                FUN=function(vec){mean(vec,na.rm=T)})

#p=ggplot()+geom_boxplot(data=tmp,aes(x=len,
#  fill=intens,y=prop),position="dodge")

tmp2 = aggregate(list(prop=tmp$prop),
                 by=list(len=tmp$len,
                         intens=tmp$intens),
                 FUN=function(vec){mean(vec,na.rm=T)})
tmp2$length = sapply(strsplit(x=as.character(tmp2$len),split = "en_"),function(el){el[2]})                
tmp2$length = factor(tmp2$length,
                     levels=c('Inf',as.character(1:(groupFrom-1)),paste(groupFrom,'andMore',sep="") ))
#tmp2$log = log10(tmp2$prop)
p=ggplot()+geom_col(data=tmp2,
                    aes(x=length,
                        y=prop,
                        fill=intens,
                        group=intens),
                    position="dodge")+
  xlab('Shortest-Path length from basal to top species')+
  ylab('Average proportion across networks')+
  scale_fill_discrete(name="Land use intensity")+
  theme_bw()+
  theme(text = element_text(size=23))

setwd(saveDir)
png('Histogram_SP_lengths.png',width=1300,height=1000)
print(p)
dev.off()

