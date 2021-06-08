#######################
## Priorisation of ha
#######################
source("code/functions_OI.R",echo=F)
pacman::p_load(sf,tidyverse,data.table,raster) 


## Choose guilds
# GUILDS <- paste0("G",c(5:8,101))  # trame humide prio
# GUILDS <- paste0("G",c(14:16,102)) # trame seche prio
# GUILDS <- paste0("G",c(2,3,4,101)) # guildes non-prio humides
GUILDS <- paste0("G",c(10,12,13,22,102)) # guildes non-prio sèches
# GUILDS <- paste0("G",c(17,19,101,102)) # trame forêts
# GUILDS <- paste0("G",c(20,101,102)) # trame prairies seches altitude


## Creat matrix of overlap among guilds
LOCATION <- "C://Dossier_Ervan/R/"

# Receiving data sets
grid_sf <- LOAD("grid100","data")
PRE <-GRID <- setDT(grid_sf)[,c("grid.id","CNHA","BV04")]

for (i in 1:length(GUILDS)) {
  # Check overalp among guilds
  dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0(GUILDS[i],"-"),list.files(LOCATION))])
  pred <- LOAD("qual4leaf","data",dir)
  
  # Add predictions of each guilds to receiving table (PRE)
  PRE <- merge(PRE,pred[Qp==1,c("grid.id","Qp")],by="grid.id",all.x=TRUE)  # overlap pred
  names(PRE)[grep("Qp",names(PRE))] <- GUILDS[i]
  
  GRID <- merge(GRID,pred[!is.na(ist),c("grid.id","ist")],by="grid.id",all.x=TRUE)  # overlap obs
  names(GRID)[i+3] <- paste0("ist",substr(GUILDS[i],2,4))
}


head(PRE)
head(GRID)


### Debut de la boucle par guilde
for (GLD2 in 1:length(GUILDS)) {
  
  ########################################################
  ### OVERLAP AMONG GUILDS (EXCLUSION oF TRAMES)
  ########################################################
# - Pour l'overlap interguilde : je le mesurerai dans chaque cluster de la mani?re suivante : 
# 1) on liste toutes les cat?gories d'overlap de guilde (p. ex. G5 U G6 ; G5 U G6 U G7, etc... ) dans le potentiel 
# 2) pour chacune de ces cat?gorie d'overlap, on mesure dans la qualit? observ?e le coefficient corrig? comme protocol? en f?vrier dernier
# 3) Au sein de chaque cluster, on fait une moyenne de cette overlap et on utilise cette valeur comme correction
  T1 <- Sys.time()
  GUILD <- gsub("G","",GUILDS[GLD2])
  dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0(GUILDS[GLD2],"-"),list.files(LOCATION))])
  path2data <- paste0(dir,"/data")
  
  ## Import data quality and filter to avoid ha in federal inventories
  PRED <- LOAD("qual4leaf","data",dir)
  # P <- merge(PRED,grid_sf[,c("grid.id","au","flachm","hochm","TWW","vogel")],by="grid.id",all.x=TRUE)
  P <- setDT(PRED)[quality=="pred" & au==0 & flachm==0 & hochm==0 & TWW==0 & vogel==0]

  ## Exclude TRAMES
  if (gsub("G","",GUILDS[GLD2])%in%c(101,102)) {
    P$val = NA
  } else {
  
  # # Check dataset already exists
  # if (any(grepl("prio_overlap",list.files(path2data)))) {
  #    PP <- LOAD("prio_overlap","data",dir)
  #    if (nrow(PP)!=nrow(P)) { stop("Former overlap doesn't match with PRED")
  #    P <- merge(P,PP,by="CNHA",all.x=TRUE)
  #   } else  {
  
  # 1) list all possible overalps in pred
  cols <- grep("ist",names(GRID))
  PRE[,"sum" := rowSums(.SD, na.rm=TRUE),.SDcols=cols]
  PP <- PRE[sum>0]
  summary(PRE$sum)
  PP[is.na(PP)] <- 0
  PP[,"combi":=do.call(paste,c(.SD, sep = "-")),.SDcols=cols]
  OVER.P <- unique(PP[sum>1 & get(GUILDS[GLD2])==1]$combi )  ## toutes les combi pour la guilde en question
  
  
  # 2) overlap in quality observed (IST) by overlap
  GRID[,c(cols):=lapply(.SD,function(x) {ifelse(is.na(x),NA,1)}),.SDcols=cols]
  GRID[,"sum" := rowSums(.SD, na.rm=TRUE),.SDcols=cols]
  G <- GRID[sum>0]
  G[is.na(G)] <- 0
  G[,"combi":=do.call(paste,c(.SD, sep = "-")),.SDcols=cols]
  OVER.G <- unique(G[sum>1]$combi)
  
  OVER.P[match(OVER.P,OVER.G)]
  OVER.P[is.na(match(OVER.P,OVER.G))] # 1 class in not present in observation
  
  setDT(P)[,"val":=as.numeric(0)] # val = overalp inter-guilde
  
  ## For each combination of overlap predicted, we compute the overlap in observations
  system.time(for (j in OVER.P){
    ID <- unique(PP[combi==j]$BV04) ## BV dans lesquels il y a de l'overlap PRED
    ID.cells <- unique(PP[combi==j]$grid.id)
    
    OBS  <- data.table("CNHA"=NA,"BV04"=NA,"name"=NA,"group"=NA,"guild"=NA)
    
    for (k in GUILDS[which(stringr::str_length(GUILDS)<4)]) {   ## !!!  Evite les trames (sinon overlap =1)
      ## IST  (Check overalp among guilds)
      dir2 <- paste0(getwd(),"/data/infofauna")
      FF <- list.files(dir2)[grep(paste0("_",gsub("G","",k),"_bio-idx_infofauna"),list.files(dir2))]
      IST <- fread(paste(dir2, FF,sep="/"))
      
      temp <- FILTER.OBS(gsub("G","",k))
      temp2 <- setDT(temp)[CNHA%in%IST$CNHA]  ## filter obs in IST infofauna
       
      tt <- merge(setDT(temp2)[,c("CNHA","name","group")],GRID[,c("CNHA","BV04")],by="CNHA",all.x=TRUE)
      tt <- tt[!is.na(match(tt$BV04,ID))]
      tt$guild <- k
      OBS <- rbind(OBS,tt)
    }
    OBS <- OBS[-1,]
    TAB <- dcast(OBS,BV04+name~guild,value.var="guild",fun=function(x) length(unique(x)),margins="guild")
    head(TAB)
    
    BB <- split(TAB,TAB$BV04)
    
    library(parallel)
    ncores<-5
    cl <- makeCluster(ncores)
    clusterExport(cl,varlist=c("val.poly","BB"),envir=.GlobalEnv)
    clusterEvalQ(cl,{
      library(data.table)
    })
    system.time(test<- parLapply(cl,BB, function(x) val.poly(x)))  # 2.5 sec
    stopCluster(cl)
    
    RES <- merge(GRID[,c("grid.id","BV04")],data.table(grid.id=as.numeric(names(test)),val=unlist(test)),by="grid.id",all.x=TRUE)
    RES <- within(RES,val[GRID$sum==1] <- 1)
    RES <- RES[!is.na(val),]
    summary(RES$val)
    
    ## Moyenne de l'overlap dans les bassins versants:
    MOY <- RES[,.("val"=mean(val)),by=BV04]
    summary(MOY$val)
    
    ### Assigne la valeur moyenne de l'overalp des obs par BV 
    P <- within(P,val[grid.id%in%ID.cells] <- MOY$val[match(P[grid.id%in%ID.cells,]$BV04,MOY$BV04)])
    # P[grid.id%in%ID.cells,]$val
    
    ## Lorsqu'il n'y a pas de valeurs (pas d'obs ds le BV), on assigne l'overlap au niveau national
    P[grid.id%in%ID.cells & is.na(val),]$val <- val.poly(TAB)
  })  # 23 min
  
  temp <- P[,c("CNHA","val")]
  save(temp,file=paste0(dir,"/data/prio_overlap.Rdata"))
    # }# end of CHECK if overlap already exists
  } # end of TRAME exlcusion
    
  # # Check overalp within hectares
  # GRID[,"sum":= rowSums(GRID[,2:ncol(GRID)],na.rm=TRUE)]
  # ## Nbr hectares unique entre guildes
  # nrow(GRID[GRID$sum>0,])  # 9145 ha
  
  # # Check overalp within hectares
  # PRE$sum <- NULL
  # PRE[,"sum":= rowSums(PRE[,3:ncol(PRE)],na.rm=TRUE)]
  # ## Nbr hectares unique entre guildes
  # nrow(PRE[PRE$sum>0,])  # 108264 ha
  # 
  # P <- merge(P,setDT(PRE)[sum>0,c("grid.id","sum")],by="grid.id",all.x=TRUE)
  # # table(is.na(P$sum))
  cat(paste0("GUILD ",GUILD," - Step 1 - overlap: done in ",round(difftime(Sys.time(),T1,units="mins"),2)," min."), "\n")
   
  ########################################################
  ### Least-cost path instead of euclidean distance to NN
  ########################################################
  pacman::p_load(sf,sp,FNN,tidyverse,data.table,foreach,doParallel,gdistance) 

  # ## Check dataset already exists
  if (any(grepl("prio_connect",list.files(paste0(dir,"/data"))))) {
   PP <- LOAD("prio_connect","data",dir)
   # load(paste0(dir,"/data/prio_connect.Rdata"))
   if (nrow(PP)!=nrow(P)) { stop("Connectivity doesn't match with PRED")}
   P1 <- merge(P,PP,by="CNHA",all.x=TRUE)
   } else  {
  
  ## Read Polygons
  dir2 <- "C://Dossier_Ervan/R/INFOFAUNA"
  if (GUILD %in% c(5:9,14:16)) {
  FF <- list.files(dir2)[grep(paste0("_",GUILD,"_cHull"),list.files(dir2))]
  } else {
  FF <- list.files(dir2)[grep(paste0("_",GUILD,"_cHull.shp"),list.files(dir2))]  
  }
  pol <-st_read(paste(dir2, FF,sep="/"))
  pol$area <- as.numeric(st_area(pol))
  st_geometry(pol) <- "geometry"
  pol$npts <- mapview::npts(pol, by_feature = TRUE)
  summary(pol$npts)
  
  pol2 <- st_cast(pol,"LINESTRING")
  system.time(pts <- st_line_sample(pol2,density=1/200,type="regular")) ## 3sec
  table(st_is_empty(pts))
  MPTS <- st_cast(pol,"MULTIPOINT")
  N1 <- mapview::npts(MPTS)
 
  if (N1>1000){
  MPTS$geometry[which(pol$npts>=10)] <- pts[which(pol$npts>=10)]
  N2 <- mapview::npts(MPTS)
   if(N1<N2) { 
    cat("More points after simplification")
    stop()}
  }  
  st_geometry(MPTS) <- "geometry"
  # table(st_is_empty(MPTS))
  # head(MPTS[st_is_empty(MPTS)==TRUE,])
  # summary(MPTS[st_is_empty(MPTS),]$area)
  
  P1 <- P
  st_geometry(P1) <- "centro"
  # P1$centro <- st_geometry(st_centroid(P1))
  if(any(st_is_empty(P1))){
  P1 <- P1[!st_is_empty(P1),]
  }
  
  ### Clusterisation des ha qualites -> creation d'un cluster tous les 500-750 m
  dist.nn <- knnx.dist(st_coordinates(MPTS)[,1:2],query = st_coordinates(st_as_sf(P1))[,1:2],algo = 'kd_tree',k =1) # Ne garder que les points a max 5km d'un polygone
  P1_5km <- setDT(P1)[which(dist.nn<=5000),]
  
  #### if a disaggregated dataset already exist, skip this part
  # if (any(grepl("pred_dsg1000",list.files(path2data)))) {
  #  P_d <- LOAD("pred_dsg1000","data",dir=dir) 
  #  CENTRO <- st_as_sf(P_d)
  # } else {
      # info <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/Guildes_revues_ER.xlsx",sheet=1,col_names = TRUE)
      # defrag <- info[which(info$No==GUILD),]$defrag
      st_geometry(P1_5km) <- "centro"
      P1sp <- as_Spatial(P1_5km)
      system.time(P1_dsg<-remove.duplicates(P1sp,zero=1000)) # 2.15 heures
      save(P1_dsg,file=paste0(path2data,"/pred_dsg1000.Rdata"))
      CENTRO<- st_as_sf(P1_dsg)
  # }
 
  st_crs(CENTRO) <- st_crs(MPTS)
  
  system.time(INTER <-st_intersects(st_buffer(CENTRO,5000),MPTS))  ## trouve tous les polygons ds un rayon de 5km
  
  # Enleve les centro sans polygones ds un rayon de 5 km
  A <- lapply(INTER,function(x) all(is.na(x)))
  POSI <- which(unlist(A)==TRUE)
  if (length(POSI) > 0) {
    CENTRO <- CENTRO[-POSI,]  # enl?ve les ha sans polygone ds rayon de 5 km
    INTER[POSI] <- NULL
  }
  
 #    ## Create raster friction from data (04/05/21): 1)use predictions as suitability map, 2) add observation (prob =1), 3) convert to friction: 1- pred (NA=0.99)
 grid <- raster::raster("D://SIG/data/grid100/grid100.tif")
 # grid_sf <- LOAD("grid100","data")

 # # for (G in c(2,3,4,10,12,13,17,19,20,22)){
 #  dir <- paste0("C:/Dossier_Ervan/R/",list.files("C:/Dossier_Ervan/R")[grep(paste0("G",G,"-"),list.files("C:/Dossier_Ervan/R"))])
 #  P <- LOAD("qual4leaf","data",dir)
 #  P[,"friction":=1-Wmean]  ## friction 0=suitable/1=unsuitable
 #  P <- within(P,friction[!is.na(ist)]<-0)  ## all obs have a friction of 0 (prob Wmean =1)
 #  ## Merge with Swiss extent (discard lacs)
 #  GRI <- merge(setDT(grid_sf[lac==0,c("CNHA","geometry")]),setDT(P[,c("CNHA","friction")]),by="CNHA",all.x=TRUE)
 #  GRI <-  GRI[is.na(get("friction")), "friction":=0.99] # set all grid without pred/obs to 0.99
 #  sf::st_geometry(GRI) <- "geometry"
 #  RA <- fasterize::fasterize(GRI,grid,field="friction")
 #  raster::plot(RA)
 #  raster::writeRaster(RA,file=paste0("C:/Dossier_Ervan/R/INFOFAUNA/raster_friction/","g",G,".tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
 # }
 
    ##Import raster friction Infofauna
    d <- "C:/Dossier_Ervan/R/INFOFAUNA/raster_friction"
    FF <- list.files(d)[grep(paste0("g",GUILD,".tif"),list.files(d))]
    r2 <- raster::raster(paste(d,FF,sep="/")) 
    rs <- raster::readAll(r2)
    
    st_geometry(P) <- "geometry"
    Q_ra <- fasterize::fasterize(P,grid)
      if(table(is.na(values(Q_ra)))[1]!=nrow(P)) {
        stop()
        cat("Different nbr of cells in raster & pred")
      }
    Q_ra <- projectRaster(Q_ra,crs="+init=epsg:2056",method="ngb")
    writeRaster(Q_ra,paste0("D://SIG/SIG_BAFU/Produits_IE/projection/G",GUILD,"_EG_2056.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
  
  #### Compute least cost paths
  registerDoParallel(makeCluster(6)) # Change the number of core

  system.time(test <- foreach(i=1:nrow(CENTRO), .packages=c('sp','gdistance','raster','sf')) %dopar% {
    origin <- CENTRO[i,]
    target <- sf::st_cast(MPTS[INTER[[i]],],"POINT",warn=F)
    dist2poly <- tapply(as.numeric(sf::st_distance(origin,target)),target$field_1,min)
    id.poly <- which(dist2poly<=5000) # select only polygons within 5km radius
    dMAX <- ifelse(max(dist2poly[id.poly])*1.2<500,500,max(dist2poly[id.poly])*1.2) # si dist2poly < 100 => LCP ne fontionne pas.
    target <- sf::st_cast(MPTS[INTER[[i]][id.poly],],"POINT",warn=F)
    
    if(dMAX >= 8000) {NA
    }else {
      BB <- extent(st_coordinates(origin)[1]-dMAX,st_coordinates(origin)[1]+dMAX,
                   st_coordinates(origin)[2]-dMAX,st_coordinates(origin)[2]+dMAX)
      zoom <- raster::crop(rs,BB)
      zoom[is.na(zoom)]<-min(getValues(zoom),na.rm=TRUE) # removing NA so as we can get values
      tr1 <- transition(zoom, transitionFunction=mean, directions=8)
      tr1c <- geoCorrection(tr1)
      A2B <- tryCatch(A2B <- costDistance(tr1c, SpatialPoints(st_coordinates(origin)), st_coordinates(target)),error=function(e) {e})
      if (is.list(A2B)|all(A2B==Inf)) {  # if costDist returns an error or Inf (target on the other side of a lake)
        NA
      } else {
        names(A2B) <- target$field_1
        LCP <- A2B[which.min(A2B)]
        a.pol <- pol[pol$field_1%in%names(LCP),]$area  # surface des polygones IST en ha
        
        #sum(LCP*LCP/sum(LCP))/sum(a.pol*a.pol/sum(a.pol))  # formule isolation InfoFauna FOIREUX POUR LA CONNECTIVITY
        dist.weight<-dnorm(seq(0,max(LCP),by=100),sd=1250)/max(dnorm(seq(0,max(LCP),by=100),sd=1250)) # Kernel -> sd = 2*500m (en version cost =~2*625m)
        if(sum(LCP)==0){
          sum(pol$area)  # si le point est inclus ds le polygon (LCP=0), renvoie 0
        } else {
          sum(a.pol*dist.weight[round(LCP/100)])
        }
      }}
  }  )  # 18 min 
  
    
  P1_dsg <- CENTRO 
  P1_dsg$connectivity<-unlist(test)
  summary(P1_dsg$connectivity)
  P1_dsg$connectivityOriginal<-P1_dsg$connectivity
  P1_dsg$connectivity<-P1_dsg$connectivity^0.25 # transformation to keep information at a 0-1 scale
  P1_dsg <- within(P1_dsg,connectivity[connectivity> round(quantile(connectivity,na.rm=TRUE,0.95),2)] <- round(quantile(connectivity,na.rm=TRUE,0.95),2)) # remplace ouliers par le 95e percentile !!! je serais plus conservatif et garderait le 70 percentile -> ce qui est loin n'est pas informatif
# 
  # PUP <- cbind(as.data.frame(P1_dsg[,c("CNHA","connectivity")] %>% st_drop_geometry()),st_coordinates(P1_dsg))
  # PUP[is.na(PUP)] <- 1
  # write.csv(PUP, paste0("D://SIG/SIG_BAFU/priorisation/", GUILD,"_prio_dsg.csv"),row.names = F)

  # Reassign values of desagregated points (P1_dsg) to nearest points within 5km of polygons (P1_5km)
  st_geometry(P1_5km) <- "centro"
  P1_5km_index<-get.knnx(st_coordinates(P1_dsg),st_coordinates(P1_5km),k=1,algo='kd_tree')$nn.index
  P1_5km$connectivity<-NA
  P1_5km$connectivity<-P1_dsg$connectivity[P1_5km_index]
  
  P1$connectivity <- NA
  P1$connectivity[match(P1_5km$CNHA,P1$CNHA)] <- P1_5km$connectivity
  summary(P1$connectivity)
  PP <- setDT(P1)[,c("CNHA","connectivity")]
  save(PP,file=paste0(path2data,"/prio_connect.Rdata"))
  # } # loop to load "prio_connect"
  
  cat(paste0("GUILD ",GUILD," - Step 2 - least cost path: done in ",round(difftime(Sys.time(),T1,units="mins"),2)," min."), "\n")

  
  ##########################
  ##### Donnees historiques (pas filtrees par EP, mais le merge permet d'avoir des valeurs ds EG)
  ##########################
  require(data.table)
  require(sf)
  
  ## Donn?es infofauna historiques
  dir2 <- "C://Dossier_Ervan/R/INFOFAUNA"
  if (GUILD %in% c(5:9,14:16)) {
    FF <- list.files(dir2)[grep(paste0("_",GUILD,"_histo"),list.files(dir2))]
  } else {
    FF <- list.files(dir2)[grep(paste0("_",GUILD,"_bio-idx_infofauna_histo"),list.files(dir2))]  
  }
  temp <- read.csv(paste(dir2, FF,sep="/"))
  
  # add hectare ID to pictis data frame
  toKM2 <- function(X) {
    x <- substr(X,1,3) 
    y <- substr(X,5,7)
    km2 <- as.numeric(format(paste0(x,y),scientific=F))
    km2
  }
  
  data.table::setDT(P1)[,CKM2:=toKM2(CNHA),by=seq_len(nrow(P1))]
  P <- merge(setDT(P1),temp[,c("CKM2","BIOIDX_TXG")],by="CKM2",all.x=TRUE)
  cat(paste0("GUILD ",GUILD," - Step 3 : done in ",round(difftime(Sys.time(),T1,units="mins"),2)," min."), "\n")
  
  #########################################################################
  ### Compute final index
  setDT(P)[,"guild_overlap":=val/max(val)]
  P[,"connect":=connectivity/max(connectivity,na.rm=TRUE)]
  
  # Lissage par moving window
  # prio2 <- merge(setDT(P),setDT(grid_sf)[,c("CNHA","grid.id","geometry")],by="CNHA",all.x=TRUE)
  grid <- raster::raster("D://SIG/data/grid100/grid100.tif")
  st_geometry(P) <- "geometry"
  r <- fasterize::fasterize(P,grid,field="connect") # quality
  
  MW <- matrix(1,11,11)
  MF <-getValues(focal(r, w=MW, mean,na.rm=TRUE))
  df <- data.table("grid.id"=getValues(grid),"connectivity"=MF)
  P <- merge(setDT(P)[,-c("geometry","connectivity")],setDT(df),by="grid.id",all.x=TRUE)
  
  
  P[,"historic_quality":=BIOIDX_TXG/max(BIOIDX_TXG,na.rm=TRUE)]
  
  P[,"env_suitability":=Wmean-min(Wmean)]
  P[,"env_suitability":=env_suitability/max(env_suitability)]
  
  mysum <- function(x){sum(x, na.rm=TRUE)}
  P[,"consensus":= rowSums(.SD,na.rm=TRUE),.SDcols=c("guild_overlap","connectivity","historic_quality","env_suitability")]
  summary(P$consensus)
  table(is.na(P$consensus))
  
  cols <- c("guild_overlap","connectivity","historic_quality","env_suitability","consensus")
  P <- setDT(P)[,(cols) :=round(.SD,2),.SDcols=cols]
  P4 <- setDT(P)[!duplicated(grid.id),c("CNHA", "BV04", "canton", "subreg","guild_overlap","connectivity","historic_quality","env_suitability","consensus","centro")]
  st_geometry(P4) <- "centro"
  save(P4,file=paste0(dir,"/data/",GUILD,"_prio.Rdata"))
  # 
  # st_geometry(P4) <- "centro"
  # PPP <- cbind(setDT(P4 %>%  st_drop_geometry()), st_coordinates(P4))
  # PPP[is.na(PPP)] <- 0
  # write.csv(PPP, paste0("D://SIG/SIG_BAFU/priorisation/", GUILD,"_prio.csv"),row.names = F)
  #  
  # P4 <- st_transform(P4,2056)
  # P4 <- cbind(as.data.frame(P4 %>% st_drop_geometry()), st_coordinates(P4))
  # P4[is.na(P4)] <- 0
  # write.csv(P4, paste0("D://SIG/SIG_BAFU/priorisation/", GUILD,"_prio_2056.csv"),row.names = F)
   
   
  cat(paste0("GUILD ",GUILD," done in ",round(difftime(Sys.time(),T1,units="mins"),2)," min."),"\n") 
}



