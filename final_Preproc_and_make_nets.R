# ADAPT directories 
repoDir = "C:/Users/user/pCloud local/boulot/Github/tetrapods_trophicNet_VS_anthropization/"
setwd(repoDir)
source('functions_to_source.R')
occDir = "C:/Users/user/pCloud local/boulot/data/GBIF IUCN tetrapods/"
saveDir = paste(occDir,'resultsRepro/',sep="")

load(file="raw_data")

#####
# Create raster of cell ids
#####

rCells = rLandUse
values = getValues(rCells)
rCells[values%in%c(11,12,13,90)] = NA
toNum = !is.na(rCells[])
rCells[toNum] = as.integer(1:sum(toNum))

#####
# Make study cells data.frame
#####

### Constraints on kept cells
# The first decile of sampling effort values across presence cells is 
# taken as a minimum threshold to declare species absence when it is not observed
minEffortQuantileForAbsence = 0.1
# Ratio of all GBIF species having a certain local presence OR absence status
taxoCover = .7
# Minimum observed richness
minRich = 20

OCC = merge(OCC,spTab[,c('file','Class')],by="file",all.x=T)

samp = aggregate(list(samplingEffort=rep(1,dim(OCC)[1])),
                 by=list(cell=OCC$cell,class=OCC$Class),sum)

nSp = dim(spTab)[1]
allCells = rCells[]
cells = as.data.frame(rasterToPoints(rCells))
colnames(cells)[3]="cell"
cells$region = raster::extract(shpWgs_raster,cells[,c('x','y')])
cells$landUse = raster::extract(rLandUse,cells[,c('x','y')])
spTab$nAppear = 0
cellsVsSpecies = Matrix(0,dim(cells)[1],nSp,sparse=T)
rownames(cellsVsSpecies) = cells$cell;colnames(cellsVsSpecies) = spTab$name
cdiucn = ((OCC$hasIUCNRange & OCC$Is_In_IUCN)  | (!OCC$hasIUCNRange))
cells$nCertainSpecies = nSp
for(i in 1:nSp){
  occ = OCC[cdiucn & OCC$file==spTab$file[i],,drop=F]
  # sampling efforts of presence cells
  effortsWherePresent=samp$samplingEffort[samp$class==class & samp$cell%in%unique(occ$cell)]
  spTab$effortThreshold[i] = quantile(effortsWherePresent,
          probs=minEffortQuantileForAbsence)
  PresOrAbs = union(unique(occ$cell) , 
      samp$cell[samp$class==class & samp$samplingEffort>spTab$effortThreshold[i]])
  iucnAbsenceCells = iucnAreasLists[spTab$Class[i]][[1]][spTab$file[i]][[1]]
  if(!is.null(iucnAbsenceCells) & !is.na(iucnAbsenceCells[1])){
    uncertainCells = setdiff( iucnAbsenceCells, PresOrAbs )
  }else{
    uncertainCells = setdiff( allCells , PresOrAbs )}
  spTab$nUncertainCells[i] = length(uncertainCells)
  cd1 = cells$cell%in%uncertainCells
  cells$nCertainSpecies[cd1] =  cells$nCertainSpecies[cd1] - 1
  
  # Site X species Matrix (pres/abs)
  cd = rownames(cellsVsSpecies)%in%unique(occ$cell)
  spTab$nAppear[i] = sum(cd)
  if(sum(cd)>0){cellsVsSpecies[cd,spTab$name[i]] = 1}
  
  if(i/5==round(i/5)){
    flush.console()
    cat('\r Process...',round(1000*i/nSp)/10,'%')
  }
}
cells = cells[!is.na(cells$landUse) & !is.na(cells$region),]
# remove BLACK and STEPPIC as they lack 
# important levels for our design  
cells = cells[cells$region%in%c(1,3,4,6,7,9) & !cells$landUse%in%c(71,72,74,75,80,90),]
cellsVsSpecies = cellsVsSpecies[as.character(cells$cell),]
cells$richness = rowSums(cellsVsSpecies)
cells = cells[cells$nCertainSpecies>(taxoCover*nSp) & cells$richness>minRich,]
cellsVsSpecies = cellsVsSpecies[as.character(cells$cell),]
cells$landGroup = paste(cells$region,'_',cells$landUse,sep="")
cells = merge(cells,keptLandUses[,c('code','intens','landSys')],by.x='landUse',by.y='code',all.x=T)
cells$regAndLandSys = paste(cells$region,'_',cells$landSys,sep="")
cells$intens = factor(cells$intens , levels=c('low','mid','high'))
cells$intens = C(cells$intens,contr.treatment)
cells$landSys = factor(cells$landSys, levels=c("forest",as.character(unique(cells$landSys)[unique(cells$landSys)!="forest"])))
cells$landSys = C(cells$landSys, contr.treatment)
cells$region = factor(cells$region , levels=c(7,unique(cells$region)[unique(cells$region)!=7]))
cells$region = C(cells$region , contr.treatment)
rownames(cells) = cells$cell

setwd(saveDir)
save(metaweb_Adj,
     metaweb_AdjAdult,
     dietNames,
     rCells,
     rLandUse,
     spTab,
     cells,
     cellsVsSpecies,
     file="preprocessed_data")

#####
# Make graphs list
#####

gList = list()
for(i in 1:dim(cellsVsSpecies)[1]){
  assemblage = colnames(cellsVsSpecies)[cellsVsSpecies[i,]>0]
  adj = metaweb_Adj[assemblage,assemblage]
  adj = t(adj) # edge goes from prey to predator
  gList[[i]]=graph_from_adjacency_matrix(adj,mode='directed')
  names(gList)[[i]]=rownames(cellsVsSpecies)[i]
  if(i/20==round(i/20)){
    flush.console()
    cat('\r Processed...',round(1000*i/dim(cellsVsSpecies)[1])/10,'%')
  }
}

setwd(saveDir)
save(gList,file = 'TrophicNetworksList')


