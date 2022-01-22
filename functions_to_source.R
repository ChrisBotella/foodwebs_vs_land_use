require(raster)
#require(CoordinateCleaner)
require(sp)
require(sf)
require(rgdal)
require(ggplot2)
require(Matrix)
require(igraph)
require(network)
require(intergraph)
require(magrittr)
require(GGally)
require(vegan)
require(npmv)
require(corrplot)
require(SDMTools)
#require(umap)
#require(ade4)


keptLandUses = data.frame(code=c(41,42,43,51,52,53,31,32,61,62,63,76,77,78,21,22,23),
                          short=c(paste(c('low','medium','high'),"intensity forest"),paste(c('low','medium','high'),"intensity grassland"),paste(c('extensive','intensive'),"perma. cropland"),paste(c('low','medium','high'),"intensity cropland"),paste(c('low','medium','high'),"intensity agri. mosaic"),paste(c('low','medium','high'),"intensity settlement")),
                          color=c('khaki2','khaki3','khaki4','chartreuse','chartreuse2','chartreuse4','darkorchid1','darkorchid4','gold1','gold3','gold4','tomato1','tomato3','tomato4','gray75','gray50','grey35'),
                          intens= c( 'low','mid','high','low','mid','high','low','high','low','mid','high','low','mid','high','low','mid','high' ),
                          landSys = c(rep('forest',3),rep('grassland',3),rep('perm. cropland',2),rep('cropland',3),rep('agri. mosaic',3),rep('settlement',3)))
regions = data.frame(region=1:12,name=c('ALPINE','ANATOLIAN','ARTIC','ATLANTIC','BLACK SEA','BOREAL','CONTINENTAL','MACARONESIA','MEDITERRANEAN','OUTSIDE','PANNONIAN','STEPPIC'))

ebdCompo = list(c("pApexMeta"),c("pBasalMeta","pBasal"),c("dirCon"), 
                c("omniProp","omniLev"),c("maxPath","meanPath","sdPath"),c("modul","meanShortDist"),
                c('patchAntiArea','proxToBorder','divLandUse'))
names(ebdCompo) = c('apex','basal','connectance','omnivory','vertiChains','compart','fragment')
classes = c('Amphibia','Reptilia','Mammalia','Aves')

#####
# Plot
#####


panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE,breaks="fd")
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}


multiplot <- function(plots=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  numPlots = length(plots)
  print(numPlots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#####
# Foodweb metrics functions
#####

clusteringCoef = function(adjmat){
  tmp = data.frame(n=rep(0,2*dim(adjmat)[1]),ones=0)
  for(k in 1:dim(adjmat)[1]){
    ids = which(adjmat[k,]>0)
    if(length(ids)>1){
      pairs = expand.grid(p1=ids,p2=ids)
      pairs = pairs[pairs$p1<pairs$p2,,drop=F]
      vTmp = sapply(1:dim(pairs)[1],function(j){adjmat[pairs$p1[j],pairs$p2[j]]})
      tmp$ones[k] = sum(vTmp)
      tmp$n[k] = length(vTmp)
    }
  }
  for(l in 1:dim(adjmat)[1]){
    ids = which(adjmat[,l]>0)
    if(length(ids)>1){
      pairs = expand.grid(p1=ids,p2=ids)
      pairs = pairs[pairs$p1<pairs$p2,,drop=F]
      vTmp = sapply(1:dim(pairs)[1],function(j)adjmat[pairs$p1[j],pairs$p2[j]])
      tmp$ones[k+l] = sum(vTmp)
      tmp$n[k+l] = length(vTmp)
    }
  }
  return(sum(tmp$ones)/sum(tmp$n))
}

short.weighted.trophic.levels = function(g){
  # The short weighted trophic level (SWTL) 
  # is the average (directed) distance to the basal species (who have no prey)
  
  # NB: it is computed only over the basal species that
  # can be reached from the focal species following
  # edges directions
  
  # Generalisation to non-connected graphs:
  # the SWTL of a species i
  # is the average of the distances between i 
  # and the basals species that are connected with i
  # because g is not necessarily connected
  SWtrophicLevels = NULL
  
  membership=igraph::components(g)$membership
  for(m in unique(membership)){# For each connected subgraph
    originalIds = which(membership==m)
    sub_g= induced_subgraph(g,vids = originalIds)# get the subgraph
    
    Adj = as_adj(sub_g)
    
    # ids of the basal species (who have no prey)
    basals = which(rowSums(Adj)==0)
    # Compute shortest distances between species 
    D = igraph::distances(sub_g)
    for(i in 1:dim(Adj)[1]){
      # Short Weighted Trophic Level of species i
      SWtrophicLevels[originalIds[i]] = mean(D[i,basals])
    }
  }
  return(SWtrophicLevels)
}

short.weighted.trophic.levels = function(g){
  # The short weighted trophic level (SWTL) 
  # is the average (directed) distance to the basal species (who have no prey)
  
  # /!\ By convention 
  # i predates j iif as_adjacency_matrix(g)[i,j]=1
  
  # NB: The SWTL of a species is computed only over 
  # basal species that can be reached from it 
  # following edges directions. 
  # It can be thus computed on non-connected
  # graphs
  
  SWtrophicLevels = NULL
  
  Adj = as.matrix(as_adj(g))
  # Remove Cannibalism to integrate basal species that are cannibals
  diag(Adj)=0
  
  # ids of the basal species (who have no prey)
  basals = which(rowSums(Adj)==0)
  # Compute shortest distances between species 
  D = igraph::distances(g,mode='out')
  
  if(length(basals)>0){
    for(i in 1:dim(Adj)[1]){
      # Short Weighted Trophic Level of species i
      # toDisp=which(is.finite(D[i,]))
      #D[i,toDisp]
      #Ord2 = which(is.finite(D[toDisp[toDisp!=i],]))
      DirectedDistToBasals = D[i,basals]
      DirectedDistToBasals = DirectedDistToBasals[is.finite(DirectedDistToBasals)]
      SWtrophicLevels[i] = mean(DirectedDistToBasals)
    }
    return(SWtrophicLevels)
  }else{
    print('Short Weighted Trophic Level can not be computed, there is no basal species')
    return(NULL)
  }
}

omnivory.levels = function(g){
  # Standard Deviation of the Trophic Level (Kefi) of the preys
  
  # By convention it is 0 for basal species
  # and species with only one prey.
  
  trophicLevels = trophiclevel(g)
  Adj =  as.matrix(as_adj(g))
  omnivoryLevels = NULL
  for(i in 1:dim(Adj)[1]){
    # Omnivory level of each species
    if(sum(Adj[i,]>0)>1){
      omnivoryLevels[i] = sd(trophicLevels[Adj[i,]>0])
    }else{
      omnivoryLevels[i] = 0 # For basals and those with only one prey
    }
  }
  return(omnivoryLevels)
}

mean.shortest.path.length = function(g){
  # Mean shortest directed path length between all oriented pairs of nodes 
  # discarding non-connected oriented pairs and self distance
  
  D = igraph::distances(g,mode="out")
  diag(D)=Inf
  tmp = as.vector(D)
  return(mean(tmp[is.finite(tmp)]))
}

basal.species.proportion=function(g){
  Adj = as.matrix(as_adj(g))
  basals = which(rowSums(Adj)==0)
  return(length(basals)/dim(Adj)[1])
}

intermediary.species.proportion=function(g){
  Adj = as.matrix(as_adj(g))
  inter = which(rowSums(Adj)>0 & colSums(Adj)>0)
  return(length(inter)/dim(Adj)[1])
}

top.species.proportion=function(g){
  Adj = as.matrix(as_adj(g))
  top = which(colSums(Adj)==0)
  return(length(top)/dim(Adj)[1])
}

## From Sonia Kefi's code
## C++ to R translation
## Warning: works only when basal species exist (i.e. some nodes with out-degree=0)
### NOUVEAU CODE TL corrige par Christophe
# Ce code met TL = 1 aux noeuds qui n'ont pas de proies, y compris les isol?s
trophiclevel <- function(g,nIter=10){
  adjmat <- as.matrix(as_adj(g))
  N <- nrow(adjmat)
  TLvec <- rep(1,N)
  for(rep in 1:nIter){
    TLtemp <- rep(0,N)
    for(i in 1:N){
      temp1 <- sum(adjmat[i,])  ## temp1 contains the number of prey of species i
      if(temp1>0){ ## If species i has at least one prey
        ## calculate the mean TL of the preys of species i
        TLtemp[i] = sum( adjmat[i,] * TLvec ) / temp1
      }
    }
    TLvec <- 1 + TLtemp ## TL = 1 + av TL of the prey
  }
  TLvec = TLvec - min(TLvec) +1
  return(TLvec)
}

### NOUVEAU CODE SPARSE MARC
trophicLevelMackay <- function(G) {
  A = as.matrix(get.adjacency(G))
  names_loc = rownames(A)
  u  = igraph::degree(G)
  v =  igraph::degree(G,mode ='in') -  igraph::degree(G,mode = 'out')  
  A = A[-1,-1]
  u = u[-1]
  v = v[-1]
  L = diag(u) - A - t(A)  
  L_sparse = as(L, "sparseMatrix")
  v_sparse = as(v, "sparseMatrix")  
  TL_vec = Matrix::solve(L,v)
  TL_vec = c(0,TL_vec)
  TL_vec = TL_vec - min(TL_vec)
  names(TL_vec) = names_loc
  return(TL_vec)
}

#####
# get_foodweb_metrics
#####

get_foodweb_metrics = function(gList,
                               tmpSavePath=NULL,
                               saveEvery = 10,
                               metrics=c('nV','nE','density','dirCon','modul','clust',
                                         'omniLev','omniProp','caniProp','predPerPrey','preyPerPred',
                                         'meanSWtrophLevel','meanShortestPathLength','basalProp','interProp','topProp',
                                         'vulnerabilitySD','generalitySD','sd_predPerPrey','skew_predPerPrey',
                                         'sd_preyPerPred','skew_preyPerPred','transitivity','diameter',
                                         'mean_distance','assortativity','tropLen',
                                         'mean_tropLevel','medi_tropLevel')){
  
  
  df=data.frame(nV = rep(NA,length(gList)),
                nE = NA,
                density = NA,
                dirCon = NA,
                modul = NA,
                clust = NA,
                omniLev = NA,
                omniProp = NA,
                caniProp = NA,
                predPerPrey = NA,
                sd_predPerPrey = NA,
                skew_predPerPrey = NA,
                preyPerPred = NA,
                sd_preyPerPred = NA,
                skew_preyPerPred = NA,
                meanSWtrophLevel = NA,
                meanShortestPathLength = NA,
                basalProp = NA,
                interProp = NA,
                topProp = NA,
                vulnerabilitySD = NA,
                generalitySD = NA,
                transitivity = NA,
                diameter = NA,
                mean_distance = NA,
                assortativity = NA,
                tropLen = NA,
                mean_tropLevel = NA,
                medi_tropLevel = NA,
                tropLen_MK = NA,
                mean_tropLevel_MK = NA,
                medi_tropLevel_MK = NA)
  df = df[,metrics,drop=F]
  df$name = names(gList)
  
  for(i in 1:dim(df)[1]){
    print(i)
    g = gList[[i]]
    
    n = length(V(g))
    if("nV"%in%metrics){df$nV[i] = n}
    
    if("nE"%in%metrics){df$nE[i] = length(E(g))}
    
    adjmat = t(as.matrix(as_adj(g)));rownames(adjmat)=NULL;colnames(adjmat)=NULL # Fixed transpose Adj
    
    if('density'%in%metrics){
      df$density[i] = sum( (as.vector(adjmat[lower.tri(adjmat)]) + as.vector(t(adjmat)[lower.tri(t(adjmat))])) >0 ) / (n*(n-1)/2)
    }
    
    if('dirCon'%in%metrics){df$dirCon[i] = sum(as.vector(adjmat))/(n^2)}
    
    if('modul'%in%metrics){
      adj = as.matrix(as_adjacency_matrix(g))
      adj_undir = 1. * matrix(as.vector(adj + t(adj))>0,dim(adj)[1],dim(adj)[1])
      g_undir = graph_from_adjacency_matrix(adj_undir, mode = "undirected")
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
      df$modul[i] = sum(termsMat[lower.tri(termsMat)])/(2*nLinks)
    }
    
    if('clust'%in%metrics){
      df$clust[i] = clusteringCoef(adjmat) # 0.988 correlation with Kortsch version
    }
    if(sum(c('omniLev','omniProp','caniProp')%in%metrics)>0){
      omniLevels = omnivory.levels(g)
      if('omniLev'%in%metrics){df$omniLev[i] = mean(omniLevels)}
      if('omniProp'%in%metrics){df$omniProp[i] = sum(omniLevels>0)/n}
      if('caniProp'%in%metrics){df$caniProp[i] = sum(diag(adjmat))/n}
    }
    
    # In Degrees distribution moments
    if(sum(regexpr('predPerPrey',metrics)>0)>0){
      npreds = rowSums(adjmat)
      # average count of predators
      if('predPerPrey'%in%metrics){df$predPerPrey[i] = mean(npreds)}
      if('sd_predPerPrey'%in%metrics){df$sd_predPerPrey[i] = sd(npreds)}
      if('skew_predPerPrey'%in%metrics){df$skew_predPerPrey[i] = skewness(npreds)}
    }
    
    # Out Degrees distribution moments
    if(sum(regexpr('preyPerPred',metrics)>0)>0){
      npreys = colSums(adjmat)
      # average count of preys
      if('preyPerPred'%in%metrics){df$preyPerPred[i] = mean(npreys)}
      if('sd_preyPerPred'%in%metrics){df$sd_preyPerPred[i] = sd(npreys)}
      if('skew_preyPerPred'%in%metrics){df$skew_preyPerPred[i] = skewness(npreys)}
    }
    
    # mean short weighted trophic level
    if('meanSWtrophLevel'%in%metrics){df$meanSWtrophLevel[i] = mean(short.weighted.trophic.levels(g))}
    
    # 
    if('meanShortestPathLength'%in%metrics){df$meanShortestPathLength[i] = mean.shortest.path.length(g)}
    
    # proportion of basal species
    if('basalProp'%in%metrics){df$basalProp[i] = basal.species.proportion(g)}
    # proportion of top species
    if('topProp'%in%metrics){df$topProp[i] = top.species.proportion(g)}
    # proportion of intermediary species (non basal and non top)
    if('interProp'%in%metrics){df$interProp[i] = intermediary.species.proportion(g)}
    
    if('vulnerabilitySD'%in%metrics){df$vulnerabilitySD[i] = sd(igraph::degree(g, mode = c("in")))  }# nodes in-degrees
    if('generalitySD'%in%metrics){df$generalitySD[i] = sd(igraph::degree(g, mode = c("out")))}# nodes out-degrees
    
    ## transitivity or clustering coefficient
    if('transitivity'%in%metrics){df$transitivity[i] = transitivity(g)}
    ## diameter or longuest shortest path
    if('diameter'%in%metrics){df$diameter[i]= diameter(g, directed = TRUE)}
    ## average shortest paths
    if('mean_distance'%in%metrics){df$mean_distance[i] = mean_distance(g, directed=TRUE)}
    
    ## mean distance on the undirected graph transform
    if('mean_distance_undir'%in%metrics){
      adj = as.matrix(as_adjacency_matrix(g))
      adj_undir = 1. * matrix(as.vector(adj + t(adj))>0,dim(adj)[1],dim(adj)[1])
      g_undir = graph_from_adjacency_matrix(adj_undir, mode = "undirected")
      df$mean_distance_undir = mean_distance(g, directed=TRUE)
    }
    
    
    ## assortativity or degree correlation
    if('assortativity'%in%metrics){df$assortativity[i] = assortativity.degree(g, directed=TRUE)}
    
    # Trophic levels metrics
    if(sum(c('tropLen','mean_tropLevel','medi_tropLevel')%in%metrics)>0){
      trophicLevels = trophiclevel(g)
      if('tropLen'%in%metrics){df$tropLen[i] = max(trophicLevels) - min(trophicLevels)}
      if('mean_tropLevel'%in%metrics){df$mean_tropLevel[i] = mean(trophicLevels)}
      if('medi_tropLevel'%in%metrics){df$medi_tropLevel[i] = median(trophicLevels)}
    }
    
    
    if(sum(c('tropLen_MK','mean_tropLevel_MK','medi_tropLevel_MK')%in%metrics)>0){
      trophicLevels = tryCatch(trophicLevelMackay(g),
                               error=function(cond) {
                                 message(cond)
                                 return(rep(0, length(V(g)) )  )
                               })
      if('tropLen_MK'%in%metrics){df$tropLen_MK[i] = max(trophicLevels) - min(trophicLevels)}
      if('mean_tropLevel_MK'%in%metrics){df$mean_tropLevel_MK[i] = mean(trophicLevels)}
      if('medi_tropLevel_MK'%in%metrics){df$medi_tropLevel_MK[i] = median(trophicLevels)}
    }
    
    
    if(!is.null(tmpSavePath) & i/saveEvery==round(i/saveEvery)){
      write.table(df,tmpSavePath,sep=";",row.names=F,col.names=T)
    }
  }
  write.table(df,tmpSavePath,sep=";",row.names=F,col.names=T)
  return(df)
}


#####
# Fix ssnonpartest
#####

ssnonpartest2 = function(formula, data, alpha = 0.05, test = c(0, 0, 0, 1), 
                          factors.and.variables = T) 
{
  if (!is(formula, "formula")) {
    return("Error: Please give a formula")
  }
  if (sum(test) != 1) {
    return("Error:Please specify a single test")
  }
  formula = Formula::Formula(formula)
  frame = model.frame(formula, data = data)
  if (sum(is.na(frame)) > 0) {
    return("Error: Missing Data")
  }
  groupvar.location = length(frame[1, ])
  groupvar = names(frame)[groupvar.location]
  vars = names(frame)[1:(groupvar.location - 1)]
  if (!is.factor(frame[, groupvar])) {
    frame[, groupvar] <- factor(frame[, groupvar])
  }
  levels = levels(frame[, groupvar])
  o <- order(frame[, groupvar])
  frame <- frame[o, ]
  p <- length(vars)
  a <- length(levels(frame[, groupvar]))
  exit = FALSE
  if (test[1] != 1) {
    N <- length(frame[, 1])
    ssize <- array(NA, a)
    lims <- matrix(NA, 2, a)
    for (i in 1:a) {
      ssize[i] <- length(frame[frame[, groupvar] == levels(frame[, 
                                                                 groupvar])[i], 1])
      lims[1, i] <- min(which(frame[, groupvar] == levels(frame[, 
                                                                groupvar])[i]))
      lims[2, i] <- max(which(frame[, groupvar] == levels(frame[, 
                                                                groupvar])[i]))
    }
    if (sum(ssize < 2) > 0) {
      return("Error: Each group must have sample size of at least 2")
    }
    Rmat <- matrix(NA, N, p)
    for (j in 1:p) {
      Rmat[, j] <- rank(frame[, vars[j]], ties.method = "average")
    }
    Rbars <- matrix(NA, a, p)
    for (i in 1:a) {
      for (j in 1:p) {
        Rbars[i, j] <- mean(Rmat[(lims[1, i]:lims[2, 
                                                  i]), j])
      }
    }
    Rtilda <- (1/a) * colSums(Rbars)
    Rbarovr <- (1/N) * colSums(Rmat)
    H1 <- matrix(0, p, p)
    H2 <- matrix(0, p, p)
    G1 <- matrix(0, p, p)
    G2 <- matrix(0, p, p)
    G3 <- matrix(0, p, p)
    for (i in 1:a) {
      H1 <- H1 + ssize[i] * (Rbars[i, ] - Rbarovr) %*% 
        t(Rbars[i, ] - Rbarovr)
      H2 <- H2 + (Rbars[i, ] - Rtilda) %*% t(Rbars[i, ] - 
                                               Rtilda)
      for (j in 1:ssize[i]) {
        G1 <- G1 + (((Rmat[(lims[1, i] + j - 1), (1:p)]) - 
                       Rbars[i, ]) %*% t((Rmat[(lims[1, i] + j - 1), 
                                               (1:p)]) - Rbars[i, ]))
        G2 <- G2 + (1 - (ssize[i]/N)) * (1/(ssize[i] - 
                                              1)) * (((Rmat[(lims[1, i] + j - 1), (1:p)]) - 
                                                        Rbars[i, ]) %*% t((Rmat[(lims[1, i] + j - 1), 
                                                                                (1:p)]) - Rbars[i, ]))
        G3 <- G3 + (1/(ssize[i] * (ssize[i] - 1))) * 
          (((Rmat[(lims[1, i] + j - 1), (1:p)]) - Rbars[i, 
                                                        ]) %*% t((Rmat[(lims[1, i] + j - 1), (1:p)]) - 
                                                                   Rbars[i, ]))
      }
    }
    if (det(H1) == 0 | det(H2) == 0 | det(G1) == 0 | det(G2) == 
        0 | det(G3) == 0 && test[1] != 1) {
      cat("Rank Matrix is Singular, only ANOVA test can be calculated \n")
      test = c(1, 0, 0, 0)
    }
    else {
      H1 <- (1/(a - 1)) * H1
      H2 <- (1/(a - 1)) * H2
      G1 <- (1/(N - a)) * G1
      G2 <- (1/(a - 1)) * G2
      G3 <- (1/a) * G3
      if (det(H1) == 0 | det(H2) == 0 | det(G1) == 0 | 
          det(G2) == 0 | det(G3) == 0 && test[1] != 1) {
        test = c(1, 0, 0, 0)
        cat("Rank Matrix is Singular, only ANOVA test can be calculated \n")
      }
    }
  }
  if (test[1] == 1) {
    #The ANOVA type statistic will be used in the following test \n")
  }
  if (test[2] == 1) {
    #The Lawley Hotelling type (McKeon's F approximation) statistic will be used in the following test \n")
  }
  if (test[3] == 1) {
    #The  Bartlett-Nanda-Pillai type (Muller's F approximation) statistic will be used in the following test \n")
  }
  if (test[4] == 1) {
    #The Wilks' Lambda type statistic will be used in the following test \n")
  }
  #base <- basenonpartest(frame, groupvar, vars, tests = test)
  base <- nonpartest(formula = formula,data= frame, tests = test,plots = F,permtest = F)
  if (test[1] == 1) {
    #testpval = base$pvalanova
    testpval = base$results$`P-value`[1]
  }
  if (test[2] == 1) {
    #testpval = base$pvalLH
    testpval = base$results$`P-value`[2]
  }
  if (test[3] == 1) {
    #testpval = base$pvalBNP
    testpval = base$results$`P-value`[3]
  }
  if (test[4] == 1) {
    #testpval = base$pvalWL
    testpval = base$results$`P-value`[4]
  }
  if (testpval >= alpha) {
    return(list())
  }
  if (p > a || factors.and.variables == TRUE) {
    chrichri_output = list(levels)
    if (length(levels) <= 2 && factors.and.variables == FALSE) {
      return(chrichri_output)
    }
    if (length(levels) <= 2 && factors.and.variables == TRUE) {
      exit = TRUE
    }
    if (exit == FALSE) {
      step2subsets = vector("list", a)
      for (i in 1:a){
        subsetframe <- subset(frame, frame[, groupvar] != 
                                levels[i])
        subsetframe <- droplevels(subsetframe)
        groupvarsub = names(subsetframe)[groupvar.location]
        base <- nonpartest(formula = formula,data= subsetframe, tests = test,plots = F,permtest = F)
        if (test[1] == 1) {
          #testpval = base$pvalanova
          testpval = base$results$`P-value`[1]
        }
        if (test[2] == 1) {
          #testpval = base$pvalLH
          testpval = base$results$`P-value`[2]
        }
        if (test[3] == 1) {
          #testpval = base$pvalBNP
          testpval = base$results$`P-value`[3]
        }
        if (test[4] == 1) {
          #testpval = base$pvalWL
          testpval = base$results$`P-value`[4]
        }
        
        k = a - 1
        if (a > 3) {
          if (testpval * a/k < alpha) {
            step2subsets[[i]] = levels(subsetframe[, 
                                                   groupvarsub])
          }
          else {
            step2subsets[[i]] = NA
          }
        }
        else {
          if (testpval < alpha) {
            step2subsets[[i]] = levels(subsetframe[, 
                                                   groupvarsub])
          }
          else {
            step2subsets[[i]] = NA
          }
        }
      }
      step2subsets = step2subsets[!is.na(step2subsets)]
      chrichri_output = c(chrichri_output,step2subsets)
      
      if (length(step2subsets) <= 1 && factors.and.variables == 
          FALSE) {
        return(chrichri_output)
      }
      if (length(step2subsets) <= 1 && factors.and.variables == 
          TRUE) {
        exit = TRUE
      }
    }
    if (exit == FALSE) {
      if (length(step2subsets[[1]]) <= 2 && factors.and.variables == 
          FALSE) {
        return(chrichri_output)
      }
      if (length(step2subsets[[1]]) <= 2 && factors.and.variables == 
          TRUE) {
        exit = TRUE
      }
    }
    if (exit == FALSE) {
      step2subsetcount = length(step2subsets)
      num.intersections = ((step2subsetcount - 1) * step2subsetcount)/2
      newsubsets = vector("list", num.intersections)
      k = 1
      for (i in 1:(step2subsetcount - 1)) {
        h = i + 1
        for (j in h:step2subsetcount) {
          newsubsets[[k]] = intersect(step2subsets[[i]], 
                                      step2subsets[[j]])
          k = k + 1
        }
      }
      newsubsetcount = length(newsubsets)
      sigfactorsubsets = vector("list", newsubsetcount)
      nonsigfactorsubsets = vector("list", newsubsetcount)
      for (i in 1:newsubsetcount) {
        subsetstotest = as.factor(newsubsets[[i]])
        subsetframe <- subset(frame, frame[, groupvar] %in% 
                                subsetstotest)
        subsetframe <- droplevels(subsetframe)
        groupvarsub = names(subsetframe)[groupvar.location]
        if(F){
          base <- basenonpartest(subsetframe, groupvarsub, 
                                 vars, tests = test)
          if (test[1] == 1) {
            testpval = base$pvalanova
          }
          if (test[2] == 1) {
            testpval = base$pvalLH
          }
          if (test[3] == 1) {
            testpval = base$pvalBNP
          }
          if (test[4] == 1) {
            testpval = base$pvalWL
          }
        }
        base <- nonpartest(formula = formula,data= subsetframe, 
                           tests = test,plots = F,permtest = F)
        if (test[1] == 1) {
          #testpval = base$pvalanova
          testpval = base$results$`P-value`[1]
        }
        if (test[2] == 1) {
          #testpval = base$pvalLH
          testpval = base$results$`P-value`[2]
        }
        if (test[3] == 1) {
          #testpval = base$pvalBNP
          testpval = base$results$`P-value`[3]
        }
        if (test[4] == 1) {
          #testpval = base$pvalWL
          testpval = base$results$`P-value`[4]
        }
        
        k = length(subsetstotest)
        if (a > 3) {
          if (testpval * (a/k) >= alpha) {
            nonsigfactorsubsets[[i]] = levels(subsetframe[, 
                                                          groupvarsub])
          }
          else {
            nonsigfactorsubsets[[i]] = NA
          }
          if (testpval * (a/k) < alpha) {
            sigfactorsubsets[[i]] = levels(subsetframe[, 
                                                       groupvarsub])
          }
          else {
            sigfactorsubsets[[i]] = NA
          }
        }
        else {
          if (testpval >= alpha) {
            nonsigfactorsubsets[[i]] = levels(subsetframe[, 
                                                          groupvarsub])
          }
          else {
            nonsigfactorsubsets[[i]] = NA
          }
          if (testpval < alpha) {
            sigfactorsubsets[[i]] = levels(subsetframe[, 
                                                       groupvarsub])
          }
          else {
            sigfactorsubsets[[i]] = NA
          }
        }
      }
      nonsigfactorsubsets = nonsigfactorsubsets[!is.na(nonsigfactorsubsets)]
      sigfactorsubsets = sigfactorsubsets[!is.na(sigfactorsubsets)]
      chrichri_output = c(chrichri_output,sigfactorsubsets) 
      
      if (length(sigfactorsubsets) <= 1 && factors.and.variables == 
          FALSE) {
        return(chrichri_output)
      }
      if (length(sigfactorsubsets) <= 1 && factors.and.variables == 
          TRUE) {
        exit = TRUE
      }
    }
    if (exit == FALSE) {
      if (length(sigfactorsubsets[[1]]) <= 2 && factors.and.variables == 
          FALSE) {
        return(chrichri_output)
      }
      if (length(sigfactorsubsets[[1]]) <= 2 && factors.and.variables == 
          TRUE) {
        exit = TRUE
      }
    }
    number.elements = a - 2
    for (l in 1:(a - 4)) {
      if (exit == FALSE) {
        newsubsetcount = length(sigfactorsubsets)
        rows = ((newsubsetcount - 1) * newsubsetcount)/2
        newsubsets = vector("list", rows)
        k = 1
        for (i in 1:(newsubsetcount - 1)) {
          h = i + 1
          for (j in h:newsubsetcount) {
            newsubsets[[k]] = intersect(sigfactorsubsets[[i]], 
                                        sigfactorsubsets[[j]])
            k = k + 1
          }
        }
        newsubsets = unique(newsubsets)
        if (length(nonsigfactorsubsets) > 0) {
          for (i in 1:length(nonsigfactorsubsets)) {
            for (j in 1:length(newsubsets)) {
              if (sum(!(newsubsets[[j]] %in% nonsigfactorsubsets[[i]])) == 
                  0) {
                newsubsets[[j]] = NA
              }
            }
          }
        }
        if (length(newsubsets) == 1 && is.na(newsubsets) && 
            factors.and.variables == FALSE) {
          return(chrichri_output)
        }
        if (length(newsubsets) == 1 && is.na(newsubsets) && 
            factors.and.variables == TRUE) {
          exit = TRUE
        }
        if (exit == FALSE) {
          newsubsets = unique(newsubsets)
          for (i in 1:length(newsubsets)) {
            if (length(newsubsets[[i]]) < (number.elements - 
                                           1)) {
              newsubsets[[i]] = NA
            }
          }
          newsubsets = newsubsets[!is.na(newsubsets)]
          if (length(newsubsets) == 0 && factors.and.variables == 
              FALSE) {
            return(chrichri_output)
          }
          if (length(newsubsets) == 0 && factors.and.variables == 
              TRUE) {
            exit = TRUE
          }
          if (exit == FALSE) {
            newsubsetcount = length(newsubsets)
            sigfactorsubsets = vector("list", newsubsetcount)
            nonsigfactorsubsets.new = vector("list", 
                                             newsubsetcount)
            for (i in 1:newsubsetcount) {
              subsetstotest = as.factor(newsubsets[[i]])
              subsetframe <- subset(frame, frame[, groupvar] %in% 
                                      subsetstotest)
              subsetframe <- droplevels(subsetframe)
              groupvarsub = names(subsetframe)[groupvar.location]
              if(F){
                base <- basenonpartest(subsetframe, groupvarsub, 
                                       vars, tests = test)
                if (test[1] == 1) {
                  testpval = base$pvalanova
                }
                if (test[2] == 1) {
                  testpval = base$pvalLH
                }
                if (test[3] == 1) {
                  testpval = base$pvalBNP
                }
                if (test[4] == 1) {
                  testpval = base$pvalWL
                }
              }
              base <- nonpartest(formula = formula,data= subsetframe, 
                                 tests = test,plots = F,permtest = F)
              if (test[1] == 1) {
                #testpval = base$pvalanova
                testpval = base$results$`P-value`[1]
              }
              if (test[2] == 1) {
                #testpval = base$pvalLH
                testpval = base$results$`P-value`[2]
              }
              if (test[3] == 1) {
                #testpval = base$pvalBNP
                testpval = base$results$`P-value`[3]
              }
              if (test[4] == 1) {
                #testpval = base$pvalWL
                testpval = base$results$`P-value`[4]
              }
              k = length(subsetstotest)
              if (a > 3) {
                if (testpval * (a/k) >= alpha) {
                  nonsigfactorsubsets.new[[i]] = levels(subsetframe[, 
                                                                    groupvarsub])
                }
                else {
                  nonsigfactorsubsets.new[[i]] = NA
                }
                if (testpval * (a/k) < alpha) {
                  sigfactorsubsets[[i]] = levels(subsetframe[, 
                                                             groupvarsub])
                }
                else {
                  sigfactorsubsets[[i]] = NA
                }
              }
              else {
                if (testpval >= alpha) {
                  nonsigfactorsubsets.new[[i]] = levels(subsetframe[, 
                                                                    groupvarsub])
                }
                else {
                  nonsigfactorsubsets.new[[i]] = NA
                }
                if (testpval < alpha) {
                  sigfactorsubsets[[i]] = levels(subsetframe[, 
                                                             groupvarsub])
                }
                else {
                  sigfactorsubsets[[i]] = NA
                }
              }
            }
            nonsigfactorsubsets.new = nonsigfactorsubsets.new[!is.na(nonsigfactorsubsets.new)]
            sigfactorsubsets = sigfactorsubsets[!is.na(sigfactorsubsets)]
            nonsigfactorsubsets = c(nonsigfactorsubsets, 
                                    nonsigfactorsubsets.new)
            chrichri_output = c(chrichri_output,sigfactorsubsets)
            if (length(sigfactorsubsets) <= 1 && factors.and.variables == 
                FALSE) {
              return(chrichri_output)
            }
            if (length(sigfactorsubsets) <= 1 && factors.and.variables == 
                TRUE) {
              exit = TRUE
            }
          }
        }
      }
      if (exit == FALSE) {
        if (length(sigfactorsubsets[[1]]) <= 2 && factors.and.variables == 
            FALSE) {
          return(chrichri_output)
        }
        if (length(sigfactorsubsets[[1]]) <= 2 && factors.and.variables == 
            TRUE) {
          exit == TRUE
        }
      }
      number.elements = number.elements - 1
    }
  }
  return(chrichri_output)
}


#####
# Aggregate graph
#####


aggregateGraph = function(g,nodeGroups){
  edgeList = as.data.frame(as_edgelist(g, names = TRUE))
  colnames(edgeList)=c('from','to')
  grTab = data.frame(name=names(V(g)),Gr=nodeGroups)
  edgeList = merge(edgeList,grTab,by.x="from",by.y="name",all.x=T)
  colnames(edgeList)[colnames(edgeList)=="Gr"] = 'GrFrom'
  edgeList = merge(edgeList,grTab,by.x="to",by.y="name",all.x=T)
  colnames(edgeList)[colnames(edgeList)=="Gr"] = 'GrTo'
  grEdgeList = aggregate(list(count=rep(1,dim(edgeList)[1])),
                         by=list(from=edgeList[,'GrFrom'],
                                 to=edgeList[,'GrTo']),
                         sum)
  grEdgeList = grEdgeList[order(grEdgeList$count),]
  aG = igraph::graph_from_edgelist(
    as.matrix(grEdgeList[,c('from','to')]),directed=T)
  aG = set_edge_attr(aG,"nLinksBetweenGroups",value =grEdgeList$count)
  grCount = table(nodeGroups)
  aG = set_vertex_attr(aG,"nNodesPerGroup",value= as.numeric(grCount[names(V(aG))]))
  return(aG)
}



#####
# Convolution
#####

convolution = function(matr,filter){
  
  dimR = dim(filter)[1];dimC = dim(filter)[2]
  centR = floor(dimR/2)+1;centC = floor(dimC/2)+1
  
  res = matr
  res[] = 0
  for(i in 1:dim(matr)[1]){
    for(j in 1:dim(matr)[2]){
      begR = max(1,centR-i+1)
      endR = min(dimR,centR+dim(matr)[1]-i)
      begC = max(1,centC-j+1)
      endC = min(dimC,centC+dim(matr)[2]-j)
      subFilt = as.vector(filter[begR:endR,begC:endC])
      
      bacAddR = centR-begR
      forAddR = endR-centR
      bacAddC = centC-begC
      forAddC = endC-centC
      
      res[i,j] = sum(as.vector(matr[(i-bacAddR):(i+forAddR),(j-bacAddC):(j+forAddC)])* subFilt)
    }
  }
  return(res)
}
