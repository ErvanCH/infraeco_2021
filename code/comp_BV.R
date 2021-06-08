pacman::p_load(rgdal,raster,ade4,adehabitatHR,adehabitatMA,ecospat,sf,data.table,tidyverse)
source("code/functions_OI.R")

dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",GUILD,"-"),list.files(LOCATION))])

## load
env.pot <- LOAD("envpot","data",dir)
env.pot$bioregion <- env.pot$subreg
BIOR <- split(env.pot,env.pot$bioregion)

## Load enviro variable
TT2 <- readxl::read_xlsx("data/var_enviro_selected.xlsx")
VAR <- TT2$label[!is.na(TT2[,stringr::str_which(names(TT2),paste0("G",GUILD,"$"))])]

AXES <-data.frame(a=c(1,1,2),b=c(2,3,3))
R = 100
mask <- ascgen(SpatialPoints(cbind((0:(R))/R, (0:(R)/R))),
               nrcol = R - 2, count = FALSE)

# source('~/R/functions_ER.R')
# system.time(L <- lapply(BIOR, function(x) pca.rasterize(x,axes=AXES,mask=MASK,R=100,z.th=0.01,save_raster=T))) # 8 min

dest.folder=paste0(dir,"/rasterBV/")

library(parallel)
ncores<-detectCores()-1
cl <- makeCluster(ncores,outfiles = 'dest.folder' )

clusterExport(cl,varlist=c('pca.rasterize',
                           'rasterizeBV',
                           'env.overlap',
                           'AXES',
                           'mask',
                           'R'),
              envir=.GlobalEnv)

clusterEvalQ(cl,{
  library(rgdal)
  library(raster)
  library(ade4)
  library(adehabitatHR)
  library(adehabitatMA)
  library(ecospat)
  library(sf)
  library(data.table)
  library(tidyverse)
  library(CircStats)
})

system.time(L <- parLapply(cl,BIOR, pca.rasterize,VAR=VAR,axes=AXES,mask=mask,R=R,z.th=0.01,save_raster=FALSE)) #4 min
save(L,file=paste0(dest.folder,GUILD,"_comp_bv.RData"))

## Clustering par bioregion
# Clustering
L <- LOAD("comp_bv","rasterBV",dir)
CENTERS <- min(sapply(L,function(x) length(unique(x$BV1))))
test2 <- lapply(L,CLUST,CENTERS)
test <- rbindlist(test2)
test$BASIS_NR <- as.numeric(substr(test$BV04,2,8))
test$BV04 <- as.numeric(substr(test$BV04,2,8))
test$bioreg <- rep(names(test2),as.numeric(sapply(test2,nrow)))

## Ensure that clusters don't have single BV (N>5)
load("D://SIG/data/Bassins_versants/BV_unifies_07.20.Rdata")  ## data = TT
BV <- dplyr::right_join(test,TT[,c("BV04","subreg")])  ## add BV without clusters
st_geometry(BV) <- "geometry"
BV$area <- as.numeric(sf::st_area(BV))/10000
BV$CLUST <- as.factor(paste(BV$subreg,BV$clust,sep="_"))

BV <- within(BV,CLUST[is.na(clust)]<-NA)
table(BV$CLUST)


## Check isolated BV (not surrounded by any other member of the same cluster)
# DEBUG: 
# GUILD = 2; i=20751
BV.ID <- unique(PRED$BV04)
st_geometry(BV) <- "geometry"
BV$outlier <- NA
BV <- within(BV, outlier[BV04%in%BV.ID] <- 0)
    for (i in BV.ID){
      focal <- BV[BV$BV04==i,]
      SUBREG <- unique(BV[BV$BV04==i,]$subreg)
      BVinSUBREG <- BV[BV$subreg==SUBREG & BV$BV04%in%BV.ID & !is.na(BV$CLUST),]
      idx <- unlist(st_intersects(focal,BVinSUBREG)) # look for direct neighbors in the same subregion
      idx <- idx[-which(idx==which(BVinSUBREG$BV04==i))] # remove focal BV from list
      ttt <-  match(BV[which(BV$BV04==i),]$CLUST,BVinSUBREG$CLUST[idx]) # check if there is at least 1 BV from the same cluster
      if (is.na(ttt)) {
        if (all(is.na(BVinSUBREG$CLUST[idx]))) { # if focal BV is surrounded by BVs with NA CLUST, looks for 6 nearest neighbors
          K <- ifelse(nrow(BVinSUBREG)<6,nrow(BVinSUBREG),6)
          T <- suppressMessages(nngeo::st_nn(focal,BVinSUBREG,k=K,progress = FALSE))
          BV <- within(BV,CLUST[which(BV$BV04==i)] <- names(which.max(table(BVinSUBREG[unlist(T),]$CLUST))))
          BV$outlier[which(BV$BV04==i)] <- 1
        } else {
          B <- BVinSUBREG[idx,] 
          if (length(B$CLUST)>0) {
           BV <- within(BV,CLUST[which(BV$BV04==i)] <- names(which.max(table(B$CLUST))))
           BV$outlier[which(BV$BV04==i)]<- 1
           }
        }
    } # end of loop

    }
BV <- within(BV,CLUST[is.na(clust) & !BV04%in%BV.ID]<-NA)
table(BV$outlier)

# ### Check result
# library(RColorBrewer)
# st_geometry(BV) <- "geometry"
# BV2 <- aggregate(BV,by=list(BV$CLUST),FUN=unique,do_union=F)
# sf::st_geometry(BV2) <- "geometry"
# Ncol <- as.numeric(tapply(BV$CLUST,factor(BV$bioreg),function(x) length(unique(x))))
# 
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# COL1=sample(col_vector, sum(Ncol))
# CLUSTERS <- BV2[,1]
# names(CLUSTERS) <- c("Cluster ID","geometry")
# CLUST <- mapview::mapview(CLUSTERS,zcol="Cluster ID",col.region=COL1,legend=F,alpha.region=0.8,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(CLUSTERS,feature.id=F),homebutton = FALSE)
# CLUST + mapview::mapview(BV,alpha.region=0.2)


## Nbr de BV par cluster
COUNT <- setDT(BV)[,.N[],by=CLUST]
summary(COUNT$V1)
table(COUNT$V1<10)  # 10 cluster avec moins de 5 BV 
COUNT$V1<10

## Le calcul ne doit se faire que sur les BV présents ds PRED
BV$col <- as.numeric(0)
BV <- within(BV,col[BV04%in%BV.ID & !is.na(clust)] <- 1)
BV <- within(BV,col[(BV04%in%BV.ID & CLUST%in%COUNT$CLUST[COUNT$V1<5]) | (BV04%in%BV.ID & is.na(clust))]<-2)
table(BV$col)
st_geometry(BV) <- "geometry"
BV$centro <- st_geometry(st_centroid(BV))
st_geometry(BV) <- "centro"
SINGLEBV <- na.omit(unique(BV[BV$col==2,]$BV04))

if(length(SINGLEBV)>0) {
  for (i in 1:length(SINGLEBV)) {
    SUBREG <- unique(BV[BV$BV04==SINGLEBV[i],]$subreg)
    if (nrow(BV[BV$BV04%in%BV.ID & BV$col==1 & BV$subreg==SUBREG,])>1) {
      TT <- suppressMessages(nngeo::st_nn(BV[BV$BV04==SINGLEBV[i],],BV[BV$BV04%in%BV.ID & BV$col==1 & BV$subreg==SUBREG,],k=1,progress = FALSE))
      BV <- within(BV,CLUST[BV$BV04==SINGLEBV[i]] <- BV[BV$BV04%in%BV.ID & BV$col==1 & BV$subreg==SUBREG,]$CLUST[as.numeric(TT)])
      }
  }
  table(is.na(BV$CLUST))
  BV$CLUST <- factor(BV$CLUST)
}

# # ## Check isolated BV (based on distance to centroids of cluster)
# st_geometry(BV) <- "centro"
# for ( a in unique(BV$CLUST)) {
#   centro_clust <- st_geometry(st_centroid(st_union(BV[BV$CLUST==a,])))
#   dist2centro <- unlist(nngeo::st_nn(BV[BV$CLUST==a,],centro_clust,returnDist=TRUE,progress=FALSE)[[2]])
#   d_mean <- mean(dist2centro)
#   outliers <- BV[BV$CLUST==a,]$BV04[which(dist2centro>(mean(dist2centro)+2*sd(dist2centro)))]
#   BV <- within(BV,col[BV04%in%outliers] <- 3)
#   if (length(outliers)>0){
#     for (j in outliers){
#       SUBREG <- unique(BV[BV$BV04==j,]$subreg)
#       NN <- suppressMessages(nngeo::st_nn(BV[BV$BV04==j,],BV[BV$col==1 & BV$subreg==SUBREG,],k=1,progress=FALSE))
#       BV <- within(BV,CLUST[BV$BV04==j] <- BV[BV$col==1 & BV$subreg==SUBREG,]$CLUST[as.numeric(NN)])
#       BV <- within(BV,outlier[BV$BV04==j] <- 2)
#     }}
# }
# table(BV$outlier)
# summary(tapply(BV$subreg,BV$CLUST,function(x) length(unique(x))))


# ## Check result
# library(RColorBrewer)
# st_geometry(BV) <- "geometry"
# BV2 <- aggregate(BV,by=list(BV$CLUST),FUN=unique,do_union=F)
# sf::st_geometry(BV2) <- "geometry"
# Ncol <- as.numeric(tapply(BV$clust,factor(BV$bioreg),function(x) length(unique(x))))
# 
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# COL1=sample(col_vector, sum(Ncol))
# CLUSTERS <- BV2[,1]
# names(CLUSTERS) <- c("Cluster ID","geometry")
# CLUST <- mapview::mapview(CLUSTERS,zcol="Cluster ID",col.region=COL1,legend=F,alpha.region=0.8,map.types = "Esri.WorldImagery",popup=leafpop::popupTable(CLUSTERS,feature.id=F),homebutton = FALSE)
# CLUST + BV

save(BV,file=paste0(dir,"/data/",GUILD,"_clusterBV_",format(Sys.time(), '%d-%m-%y'),".RData"))
