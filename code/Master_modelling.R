options("scipen"=100, "digits"=4)
# install.packages("pacman")
pacman::p_load(sf, data.table,tidyverse,ggplot2,reshape2,qgraph,vegan,raster,hablar)
source("code/functions_OI.R")
# LOCATION <- "G://R/"# location of all folders (PC Infoflora)
LOCATION <- "C://Dossier_Ervan/R/"# location of all folders (IF office)

GUILD <- 10

# Creating guild space -------------------------------------------------------


# for (GUILD in G) {
  dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",GUILD,"-"),list.files(LOCATION))])
  guild <- readxl::read_xlsx("data/Info_guildes.xlsx",sheet=1,col_names = TRUE,col_types = "guess")
  guild.info <- guild[which(guild$ID==GUILD),]
  guild.info <- guild.info %>% hablar::convert(lgl(RF,GAM,GBM,GBM_xg,ME))
  DATE <- format(Sys.time(), '%d-%m-%y')
  
  # Assign spatial info & enviro variables
  TT <- readxl::read_xlsx("data/var_enviro_by_guild.xlsx")
  VAR <- TT$label[!is.na(TT[,stringr::str_which(names(TT),paste0("G",GUILD,"$"))])]
  
  ## Load IST Infofauna & attach grid info
  load("data/grid100_sf_with_enviro.Rdata")
  FF <- list.files("data/infofauna")[grep(paste0("_",GUILD,"_bio-idx"),list.files("data/infofauna"))]
  ist <- data.table::fread(paste("data/infofauna", FF,sep="/"))
  IST <- merge(setDT(ist),setDT(grid_sf)[,c("CNHA","grid.id","geometry","BV04")],by="CNHA",all.x=T)
  st_geometry(IST) <- "geometry" 


  ### Desagregation des donnees
  if (is.na(guild.info$IST_defrag)){
    cat("!! Pas de paramètres de défragmentation !!")
  }
  system.time(IST2 <- IST.dsg(IST,guild.info$IST_defrag))


  ## Data for calibration
    ## A. Observations
  div2 <- setDT(grid_sf)[grid.id%in%IST2$grid.id,c("grid.id","CNHA","BV04","canton","subreg","geometry","centro",..VAR)]
  sf::st_geometry(div2) <- "centro"
  
  Q1 <- na.omit(data.table::setDT(div2)) 
  nrow(IST2)-nrow(Q1) # 18 ha vires
  Q1$Qobs=1

  ## B. Pseudo-absence dans tout l'espace guilde
  eg <- LOAD("EG","data",dir)
  names(eg) <- c("grid.id","geometry")
  EG <-setDT(grid_sf)[grid.id%in%c(IST$grid.id,eg$grid.id),c("grid.id","CNHA","BV04","canton","subreg","geometry","centro",..VAR)]
  EG <- EG[!duplicated(EG$grid.id),]
  NArep(EG)
  
  # summary(EG)
  Q0 <- dplyr::sample_n(na.omit(data.table::setDT(EG)[!grid.id%in%Q1$grid.id]),nrow(Q1),replace=F)
  Q0$Qobs=0
  OBS <- rbind(Q1,Q0)

  # ## Selection des variables
  # library(rgdal)
  # library(rgeos)
  # library(parallel)
  # library(raster)
  # library(perm)
  # 
  # ncores<-detectCores()-1
  # if (ncores == 0) {ncores = 1}
  # 
  # env.pres<-as.data.frame(OBS[OBS$Qobs==1,..VAR])#env data for presences
  # env.bck<-as.data.frame(OBS[OBS$Qobs==0,..VAR])#env data for background
  # to.do<-1:ncol(env.pres) #index for parallelization
  # 
  # cl <- makeCluster(ncores)
  # clusterExport(cl,varlist=c('to.do','env.bck','env.pres','var.test'),envir=.GlobalEnv)
  # clusterEvalQ(cl,{
  #   library(perm)
  # })
  # 
  # #apply non-parametric t-test between background and presences for each variable
  # system.time(VAR.test<-parLapply(cl,X=to.do,fun=var.test,Nrep=5,mcmc=1000,env.pres,env.bck))  # 10 sec
  # stopCluster(cl)
  # 
  # #results formatting
  # select.var <- data.frame(var=VAR,do.call(rbind,VAR.test))
  # write.table(select.var, file = "data/selection_var.txt",sep='\t',quote=F,row.names = F)
  # 
  # ## Generate report for variable selection
  # rmarkdown::render(input = paste0(getwd(),"/code/var_select.Rmd"), 
  #                   output_format = "html_document",
  #                   output_file = paste0(GUILD,"_var_select2.html"),
  #                   output_dir = paste0(getwd(),"/report"))
  

  ## Espace de projection
  TT2 <- readxl::read_xlsx("data/var_enviro_selected.xlsx")
  VAR <- TT2$label[!is.na(TT2[,stringr::str_which(names(TT2),paste0("G",GUILD,"$"))])]
  
  ## Projection dans l'espace-guilde (to lower computational time)
  env.pot <-setDT(grid_sf)[grid.id%in%c(IST$grid.id,eg$grid.id),c("grid.id","BV04","subreg",..VAR)]
  env.pot <- na.omit(setDT(env.pot)[!duplicated(grid.id),]) 
  
  # nrow(EG) - nrow(env.pot)  # 355 ha perdus sur la fronti?re
  # st_geometry(EG) <- "geometry"
  # mapview::mapview(EG[!EG$grid.id%in%env.pot$grid.id,],legend=F)
  
  ## Projection dans l'espace-guilde (to lower computational time)
  env.pot <- merge(env.pot,setDT(IST)[,c("grid.id","BIOIDX_TXG")],by="grid.id",all.x=T)
  env.pot$ist <- env.pot$BIOIDX_TXG
  env.pot$BIOIDX_TXG <- NULL
  
  env.pot <- merge(env.pot,setDT(OBS)[,c("grid.id","Qobs")],by="grid.id",all.x=T)
  env.pot$cal <- env.pot$Qobs
  env.pot$Qobs <- NULL
  table(duplicated(env.pot$grid.id))
  env.pot <- env.pot[!duplicated(env.pot$grid.id),]
  save(env.pot,file=paste0(dir,"/data/",GUILD,"_envpot.RData"))
  
  rm(list=setdiff(ls(), c("env.pot","OBS","VAR","GUILD","guild.info","dir","DATE","LOCATION")))
  save.image(paste0(dir,"/data/",GUILD,"_enviro4model_",DATE,".RData"))


# Modelling -------------------------------------------------------
  # dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",GUILD,"-"),list.files(LOCATION))])
  # FF <- list.files(paste(dir,"data",sep="/"))[grep("enviro4model",as.vector(list.files(paste(dir,"data",sep="/"))))]
  # IF <- file.info(paste(dir,"data", FF,sep="/"))
  # load(paste(dir,"data", FF[which.max(IF$mtime)],sep="/"))
  # guild <- readxl::read_xlsx("data/Info_guildes.xlsx",sheet=1,col_names = T,col_types = "guess")
  # guild.info <- guild[which(guild$ID==GUILD),]
  # guild.info <- guild.info %>% hablar::convert(lgl(RF,GAM,GBM,GBM_xg,ME))
  rstudioapi::restartSession(command='source("code/Modelling_part.R")')
  
  
# BV clustering   -------------------------------------------------------
  options("scipen"=100, "digits"=4)
  # install.packages("pacman")
  pacman::p_load(sf, data.table,tidyverse,ggplot2,reshape2,qgraph,vegan,raster,hablar)
  source("code/functions_OI.R")
  # LOCATION <- "G://R/"# location of all folders (PC Infoflora)
  LOCATION <- "C://Dossier_Ervan/R/"# location of all folders (IF office)
 

  ### Select guild
  G <- c(2:4,10,12,13,17,19,20,22)  # nouveaux ID
  grid_sf <- LOAD("grid100","data")
  guild <- readxl::read_xlsx("data/Info_guildes.xlsx",sheet=1,col_names = TRUE,col_types = "guess")
  
  for (GUILD in G) {
  dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",GUILD,"-"),list.files(LOCATION))])
  guild.info <- guild[which(guild$ID==GUILD),]
  guild.info <- guild.info %>% hablar::convert(lgl(RF,GAM,GBM,GBM_xg,ME))
  PRED <- LOAD("quality_predicted","data",dir)
  
  if (any(grep("clusterBV",list.files(paste0(dir,"/data"))))) {
    BV <- LOAD("clusterBV","data",dir)
  }  else {
    system.time(source("code/comp_BV.R"))
    BV<- LOAD("clusterBV","data",dir)
  }
  
  
  ## Associe les clusters aux PRED
  PRED$subreg <- NULL
  PRED <- merge(setDT(PRED),setDT(BV)[,c("BV04","CLUST","subreg","area")],by="BV04",all.x=TRUE)
  table(is.na(PRED$CLUST))
  length(unique(PRED$BV04))
  
  # Avoid predictions in National Inventories (Nat.Inv.)
  PRED <- merge(PRED,grid_sf[,c("grid.id","au","flachm","hochm","TWW","vogel")],by="grid.id",all.x=TRUE)
   
  OBS <- FILTER.OBS(GUILD)
  
# Benchmarking   -------------------------------------------------------
  guild.name <- guild[which(guild$ID==GUILD),]$name_guild_de
  th.pred <-  guild[which(guild$ID==GUILD),]$th.pred
  th.bench <- guild[which(guild$ID==GUILD),]$th.bench
  th.qprop <- guild[which(guild$ID==GUILD),]$th.qprop
  defrag <- guild[which(guild$ID==GUILD),]$defrag
  bench <- guild[which(guild$ID==GUILD),]$bench
  
  RES <- BENCH(PRED,OBS,th.pred,th.bench,th.qprop,defrag,sd.min=TRUE)
  PRED <- RES[[2]]
  res <- RES[[1]]
  if(sum(res$Qpred,na.rm=T)!=nrow(PRED[PRED$quality=="pred"])) { stop("Different number of predictions")}
  save(PRED,file=paste0(dir,"/data/",GUILD,"_qual4leaf_",format(Sys.time(), '%d-%m-%y'),".Rdata"))
  
  ### Reassign geometry
  setDT(res)[is.na(CLUST),c("Nsp_BV","bench_Nsp","Qobs", "Qpred","EB_sp_defrag", "EB_sp_defrag_weighted"):=0]
  st_geometry(res) <- "geometry"
  save(res,file=paste0(dir,"/data/",GUILD,"_ha2add_",format(Sys.time(), '%d-%m-%y'),".RData"))



# Fachsheet   -------------------------------------------------------
  
  ## Specie contribution
  # RES <- unique(data.table::data.table("group"=OBS$group,"species"=OBS$name,"bioregion"=OBS$bioregrion,"weight"=OBS$w,prop_area_ha=as.numeric(NA),prop_area_percent=as.numeric(NA)))
  # NN <- data.table::setDT(OBS)[,.(Qobs=ifelse(sum(w[!duplicated(taxonid)])>=1,1,0)),by=.(grid.id)]
  # system.time(for(i in 1:nrow(RES)){
  #   NN2 <-  data.table::setDT(OBS)[name!=RES[i,]$species,.(Qobs=ifelse(sum(w[!duplicated(taxonid)])>=1,1,0)),by=.(grid.id)]
  #   RES[i,]$prop_area_ha <- nrow(NN[Qobs==1]) - nrow(NN2[Qobs==1])
  #   RES[i,]$prop_area_percent <- (nrow(NN[Qobs==1]) - nrow(NN2[Qobs==1]))*100/nrow(NN[Qobs==1])
  # }) # 70 sec
  # names(RES) <- c("group", "speciesCODE", "weight", "prop_area_ha", "Contribution")
  # write.csv(RES,file=paste0("C://Dossier_Ervan/R/INFOFAUNA/guilde_",GUILD,"_spContribution.csv"))

  # for (GUILD in c(2:4,10,12,13,17,19,20,22)) {
  #   dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",GUILD,"-"),list.files(LOCATION))])
  #   guild <- readxl::read_xlsx("data/Info_guildes.xlsx",sheet=1,col_names = T,col_types = "guess")
  bench <- guild[which(guild$ID==GUILD),]$bench
  guild.name <- guild[which(guild$ID==GUILD),]$name_guild_de
 

  res <- LOAD("ha2add","data",dir)
  res <- col.bin(res,bench,min.size=5)
  res <- res[[1]]
   
  TT <- readxl::read_xlsx("data/var_enviro_selected.xlsx")
  VAR <- TT$label[!is.na(TT[,stringr::str_which(names(TT),paste0("G",GUILD,"$"))])]
  
  rmarkdown::render(input = "code/factsheet.Rmd",
                    output_format = "html_document",
                    output_file = paste0(GUILD,"_factsheet3.html"),
                    output_dir = "D:/SIG/SIG_BAFU/factsheets_final/")


# Fiches descriptives  ------------------------------------------------------- 
  pacman::p_load(leafgl,leaflet,leafem,leafpop,data.table,raster,sf,tidyverse)
  source("code/functions_OI.R")
  guild <- readxl::read_xlsx("data/Info_guildes.xlsx",sheet=1,col_names = TRUE)
  DATE <- format(Sys.time(), '%d-%m-%y')
  load("data/grid100_sf_with_enviro.Rdata")
  
  G <- c(5:8,14:16,101,102)
  G <- c(2:4,10,12,13,17,19,20,22)  # nouveaux ID
  
  for (GUILD in G) {
    dir <- paste0("C:/Dossier_Ervan/R/",list.files("C:/Dossier_Ervan/R")[grep(paste0("G",GUILD,"-"),list.files("C:/Dossier_Ervan/R"))])
    bench <- guild[which(guild$ID==GUILD),]$bench
    guild.name <- guild[which(guild$ID==GUILD),]$name_guild_fr
    
    ## ADD IST
    if (!any(grepl("EG_fiche",list.files(paste0(dir,"/shp"))))){
      PRED <- LOAD("qual4leaf","data",dir)
      st_geometry(PRED) <- "geometry"
      PRED <- st_transform(PRED,3857)
      grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
      Q_ra <- fasterize::fasterize(PRED,grid2,field="Qp")
      raster::writeRaster(Q_ra,paste0(dir,"/shp/",GUILD,"_EG_fiche.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
    }
    
      rmarkdown::render(input="C:/Dossier_Ervan/R/Fiches/fiche_guilde_pdf.Rmd", 
                      output_format = "word_document",
                      output_file = paste0("D:/SIG/SIG_BAFU/fiches/fiche_G",GUILD,".docx"))
  }
  
  ########################
  ## GUILDES MOBILES
  ########################
  library(leafgl)
  library(leaflet)
  library(leafem)
  library(leafpop)
  library(data.table)
  library(raster)
  library(tidyverse)
  library(sf)
  source("C://Dossier_Ervan/R/functions_ER.R")
  
  guild <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/Info_guildes.xlsx",sheet=1,col_names = TRUE)
  G <- c(25,26)
  
  for (GUILD in G) {
    dir <- paste0("C:/Dossier_Ervan/R/",list.files("C:/Dossier_Ervan/R")[grep(paste0("G",GUILD,"-"),list.files("C:/Dossier_Ervan/R"))])
    # setwd(dir)
    bench <- guild[which(guild$ID==GUILD),]$bench
    guild.name <- guild[which(guild$ID==GUILD),]$name_guild_fr
    
    ## ADD IST
    if (!any(grepl("EG",list.files(paste0(dir,"/shp"))))){
      PRED <- LOAD("qual4leaf","data",dir)
      st_geometry(PRED) <- "geometry"
      PRED <- st_transform(PRED,3857)
      grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
      Q_ra <- fasterize::fasterize(PRED,grid2,field="Qp")
      raster::writeRaster(Q_ra,paste0("shp/",GUILD,"_EG.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
    }
    
    # ## Create factsheet
    # TT <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/var_enviro_selected.xlsx")
    # VAR <- TT$label[!is.na(TT[,stringr::str_which(names(TT),paste0("G",GUILD,"$"))])]
    # rmarkdown::render(input = "C:/Dossier_Ervan/R/Fiches/pdf4christophe.Rmd",
    #                   output_format = "pdf_document",
    #                   output_file = paste0("test4CRB.pdf"))
    #                   # output_file = paste0("fiche_guilde_",GUILD,".html"))
    # 
    rmarkdown::render(input="C:/Dossier_Ervan/R/Fiches/fiche_guildes_mobiles.Rmd", 
                      output_format = "word_document",
                      output_file = paste0("C:/Dossier_Ervan/R/Fiches/Fiches/fiche_guilde_",GUILD,".docx"))
  }
  
  #### GUILD sans soll (9, 18 et 23)
  GUILD = 23
  dir <- paste0("C:/Dossier_Ervan/R/",list.files("C:/Dossier_Ervan/R")[grep(paste0("G",GUILD,"-"),list.files("C:/Dossier_Ervan/R"))])
  bench <- guild[which(guild$ID==GUILD),]$bench
  guild.name <- guild[which(guild$ID==GUILD),]$name_guild_fr
  
  ## ADD IST
  if (!any(grepl("EG_fiche",list.files(paste0(dir,"/shp"))))){
    PRED <- LOAD("qual4leaf","data",dir)
    st_geometry(PRED) <- "geometry"
    PRED <- st_transform(PRED,3857)
    grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
    Q_ra <- fasterize::fasterize(PRED,grid2,field="Qp")
    raster::writeRaster(Q_ra,paste0(dir,"/shp/",GUILD,"_EG_fiche.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
  }
  
  # ## Create factsheet
  # TT <- readxl::read_xlsx("C://Dossier_Ervan/Guildes&co/var_enviro_selected.xlsx")
  # VAR <- TT$label[!is.na(TT[,stringr::str_which(names(TT),paste0("G",GUILD,"$"))])]
  # rmarkdown::render(input = "C:/Dossier_Ervan/R/Fiches/pdf4christophe.Rmd",
  #                   output_format = "pdf_document",
  #                   output_file = paste0("test4CRB.pdf"))
  #                   # output_file = paste0("fiche_guilde_",GUILD,".html"))
  # 
  rmarkdown::render(input="C:/Dossier_Ervan/R/Fiches/fiche_guilde_sans_soll.Rmd", 
                    output_format = "word_document",
                    output_file = paste0("C:/Dossier_Ervan/R/Fiches/Fiches/fiche_guilde_",GUILD,".docx"))
  
  #### GUILD sans soll (9)
  GUILD = 9
  GUILD = 18
  dir <- paste0("C:/Dossier_Ervan/R/",list.files("C:/Dossier_Ervan/R")[grep(paste0("G",GUILD,"-"),list.files("C:/Dossier_Ervan/R"))])
  rmarkdown::render(input="C:/Dossier_Ervan/R/Fiches/fiche_guilde_9.Rmd", 
                    output_format = "word_document",
                    output_file = paste0("D:/SIG/SIG_BAFU/fiches/fiche_G",GUILD,".docx"))
  
  
  
  
  
  
  
  
  
# Leaflets finaux  -------------------------------------------------------
  
  pacman::p_load(leafgl,leaflet,leafem,leafpop,data.table,raster,sf,tidyverse)
  source("code/functions_OI.R")
  guild <- readxl::read_xlsx("data/Info_guildes.xlsx",sheet=1,col_names = TRUE)
  DATE <- format(Sys.time(), '%d-%m-%y')
  load("data/grid100_sf_with_enviro.Rdata")
  
  # # G <- c(5:8,14:16,101,102)
  G <- c(2:4,10,12,13,17,19,20,22)
    
 for (GUILD in G) {
    dir <- paste0("C:/Dossier_Ervan/R/",list.files("C:/Dossier_Ervan/R")[grep(paste0("G",GUILD,"-"),list.files("C:/Dossier_Ervan/R"))])
    bench <- guild[which(guild$ID==GUILD),]$bench
    guild.name <- guild[which(guild$ID==GUILD),]$name_guild_de
    
    res <- LOAD("ha2add","data",dir)
    st_geometry(res) <- "geometry"
    
    # Add layer to GPKG
    if (GUILD!=20) {  ## Pas de BSS pour la guilde 20
    res2 <- st_transform(res,2056)
    res2 <- col.bin(res2,bench,min.size=5)[[1]]
    res2$Gilde <- paste0("G",GUILD)
    res2 <- res2[,c("BV.id","CLUST","sp.in.BV","Nsp_BIOREG","observed_qual", "potential_qual", "Erganzungsbedarf","Gilde","geometry")]
    names(res2) <- c("WaterCatchmentID","Cluster","SpeciesInWaterCatchment","SpeciesInBiogegraphicalRegion","ObservedQuality","PotentialQuality","AdditionalSurfaceNeeded","Gilde","geometry")
    write_sf(obj = res2, dsn ="D://SIG/SIG_BAFU/GPKG/BSS.GPKG",layer=paste0("G",GUILD),layer_options = "OVERWRITE=true")
    }
    
    res2 <- st_transform(res,4326)
    V1 <- col.bin(res2,bench,min.size=5)
    res <- V1[[1]]
    st_geometry(res) <- "geometry"
    

    m = leaflet() %>%
      addProviderTiles(provider = providers$CartoDB.Positron,group="Positron",layerId="Positron") %>%
      addProviderTiles(provider = providers$Esri.WorldImagery,group="Esri",layerId="Esri") %>%
      addPolygons(data = res, group = "Erganzungsbedarf",weight=1,color="grey",opacity=0.8,fillColor = ~col,fillOpacity = 0.8,highlightOptions = highlightOptions(color = "black", weight = 2),popup = popupTable(V1[[1]][,c("BV.id","sp.in.BV","sp.in.bench","observed_qual","potential_qual","Erganzungsbedarf")],feature.id = FALSE,row.numbers = FALSE)) %>%
      addLegend(colors=V1[[2]]$col,labels=V1[[2]]$lab,group = "Erganzungsbedarf", position = "topright",opacity=0.8,title=paste0("Erganzungsbedarf [ha] (max:",max(V1[[1]]$Erganzungsbedarf,na.rm=T),")")) %>%
      addPolygons(data = res, group = "BV",weight=1,color="grey",opacity=0.8,fillColor = "transparent",fillOpacity = 0,highlightOptions = highlightOptions(color = "black", weight = 2),popup = popupTable(V1[[1]][,c("BV.id","sp.in.BV","sp.in.bench","observed_qual","potential_qual","Erganzungsbedarf")]))
    
    
    m$dependencies = c(m$dependencies,
                       mapview:::popupLayoutDependencies())
    
    ## ADD IST
    dir_if <- "C:/Dossier_Ervan/R/infraeco_2021/data/infofauna"
    FF <- list.files(dir_if)[grep(paste0("_",GUILD,"_bio-idx"),list.files(dir_if))]
    ist <- data.table::fread(paste(dir_if, FF,sep="/"))
    IST <- merge(setDT(ist),setDT(grid_sf)[,c("CNHA","grid.id","geometry","BV04")],by="CNHA",all.x=TRUE)
    st_geometry(IST) <- "geometry" 
    
    setDT(IST)[,"i2":= ifelse(BIOIDX_TXG <= median(BIOIDX_TXG,na.rm=TRUE),"Hohe Qualitat","Sehr hohe Qualitat"),by=.SD]  ## binarize BIOINDEX as per BAFU-VDC layout
    st_geometry(IST) <- "geometry"
    IST <- st_transform(IST,4326)
    col1 = colourvalues::colour_values_rgb(IST$i2,palette="greens",include_alpha = FALSE)
    col11 <- colorFactor(palette = "Greens",IST$i2)
    
    m1 <- m %>%
      addGlPolygons(data = IST, group = "Observed quality (ha)",fillColor=col1,fillOpacity=0.8) %>%
      addLegend(pal =  col11, values=IST$i2,group = "Observed quality (ha)", position = "topright",opacity=2,title="Observed quality (ha)")
    
    ## ADD SOLL
    PRED <-  LOAD("qual4leaf","data",dir)
    P <- setDT(PRED)[quality=="pred",]
    # ## Priorisation du SOLL (cf code Priorisation)
    PRIO <- LOAD("prio","data",dir)
    if(nrow(P)!=nrow(PRIO)){ stop("PRED & PRIO don't have the same length")}
    
    # Add layer to GPKG
    PRIO$Gilde <- paste0("G",GUILD)
    PRI <- merge(setDT(PRIO)[,c("CNHA", "BV04", "canton", "subreg", "env_suitability","connectivity", "historic_quality", "guild_overlap","Gilde")],setDT(P)[,c("CNHA","geometry")],by="CNHA",all.y=TRUE)
    PRI[is.na(PRI)] <- 0
    names(PRI) <- c("OBJECTID","WaterCatchmentID","Canton","BiogeographicalSubRegion","EnvironmentalSuitability","Connectivity","HistoricalQuality","GuildOverlap","Gilde","geometry")
    st_geometry(PRI) <- "geometry"
    PRI <- st_transform(PRI,2056)
    write_sf(obj = PRI, dsn ="D://SIG/SIG_BAFU/GPKG/PRIO.GPKG",layer=paste0("G",GUILD),layer_options = "OVERWRITE=true")
    
    
    P <- merge(setDT(P),setDT(PRIO)[,c("CNHA","consensus")],by="CNHA",all.x=TRUE)
    table(is.na(P$consensus))
    if (any(is.na(P$consensus))) { stop("Not all predictions have a consensus")}
    st_geometry(P) <- "geometry"
    P1 <- st_transform(P,3857)
    grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
    P_ra <- fasterize::fasterize(P1,grid2,field="consensus")
    raster::writeRaster(P_ra,paste0(dir,"/shp/",GUILD,"_PRED3857.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
    
    # col3 = colourvalues::colour_values_rgb(P$prio, palette = "inferno",include_alpha = FALSE)
    col3 <- colorNumeric(palette = "inferno",P1$consensus,na.color="#00000000",reverse=TRUE)
    
    m2 <- m1 %>%
      addRasterImage(P_ra,colors=col3,method="ngb",group = "Potential quality (consensus)",maxBytes=object.size(P_ra)) %>%
      addLegend(pal =  col3,values=P$consensus,group = "Potential quality (consensus)", position = "topright",opacity=0.8,title="Potential quality (consensus)")
    
    
    # ### Add polygon INFOFAUNA
    dir2 <- "C://Dossier_Ervan/R/INFOFAUNA"
    if (GUILD %in% c(5:9,14:16)) {
      FF <- list.files(dir2)[grep(paste0("_",GUILD,"_cHull"),list.files(dir2))]
    pol <-st_read(paste(dir2, FF,sep="/"))
    pol <- st_transform(pol,4326)
    pol$PRIO <- ifelse(pol$proxy_IFed==0,"Regional","National")
    col4 = colourvalues::colour_values_rgb(pol$PRIO,palette="rdbu",include_alpha = FALSE)
    col5 <- colorFactor(palette = "RdBu",pol$PRIO)
    m3 <- m2  %>%
      addGlPolygons(data = pol, group = "Observed quality (polygons)",popup="PRIO",fillColor=col4,fillOpacity=0.8) %>%
      addLegend(pal =  col5, values=pol$PRIO,group = "Observed quality (polygons)", position = "topright",opacity=2,title="Observed quality (polygons)")
        } else {
      FF <- list.files(dir2)[grep(paste0("_",GUILD,"_cHull.shp"),list.files(dir2))]
      pol <-st_read(paste(dir2, FF,sep="/"))
      pol <- st_transform(pol,4326) 
      
      m3 <- m2  %>%
      addGlPolygons(data = pol, group = "Observed quality (polygons)",popup="id",fillColor="green",fillOpacity=0.6) %>%
      addLegend(colors ="green",labels=" ",group = "Observed quality (polygons)", position = "topright",opacity=0.6,title="Observed quality (polygons)")
    }
    
   
    ### Finalize leaflet
    m3 <- m3 %>%
      addLayersControl(baseGroups= c("Positron","Esri"),overlayGroups = c("Erganzungsbedarf","Observed quality (ha)","Observed quality (polygons)","Potential quality (consensus)","BV"),position="topleft") %>%
      hideGroup(c("Observed quality (ha)","Potential quality (consensus)","Observed quality (polygons)"))
    # mapview::mapshot(m3, url =paste0("report/",GUILD,"_final.html"))
    mapview::mapshot(m3, url =paste0("D:/SIG/SIG_BAFU/leaflets/",GUILD,"_final.html"))
    cat(paste0("GUILD ",GUILD," done!","\n","\n"))
    # Clean up a bit the working space
 
  }
  

# # Leaflet guildes sans soll -----------------------------------------------
# 
#   for (GUILD in c(9,18,24)) {
#     GUILD=18
#     dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",GUILD,"-"),list.files(LOCATION))])
#     guild <- readxl::read_xlsx("data/Info_guildes.xlsx",sheet=1,col_names = T,col_types = "guess")
# 
#     guild.name <- guild[which(guild$ID==GUILD),]$name_guild_de
# 
#     rmarkdown::render(input = "code/factsheet_sans_SOLL.Rmd",
#                       output_format = "html_document",
#                       output_file = paste0(GUILD,"_factsheet3.html"),
#                       output_dir = "D:/SIG/SIG_BAFU/factsheets_final/")
# 
#     # Create leaflets for plausibilisation
#     pacman::p_load(leafgl,leaflet,leafem,leafpop,data.table,raster,tidyverse,sf)
# 
#     load("data/grid100_sf_with_enviro.Rdata")
#     FF <- list.files("data/infofauna")[grep(paste0("_",GUILD,"_bio-idx"),list.files("data/infofauna"))]
#     ist <- read.csv(paste("data/infofauna", FF,sep="/"))
#     IST <- merge(setDT(ist),setDT(grid_sf)[,c("CNHA","grid.id","geometry")],by="CNHA",all.x=T)
#     st_geometry(IST) <- "geometry"
#     IST <- st_transform(IST,3857)
#     grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
#     I_ra <- fasterize::fasterize(IST,grid2,field="BIOIDX_TXG")
#     col1 <- colorNumeric(palette = "viridis",IST$BIOIDX_TXG,na.color="#00000000")
# 
# 
#     m = leaflet() %>%
#       addProviderTiles(provider = providers$CartoDB.Positron,group="Positron",layerId="Positron") %>%
#       addProviderTiles(provider = providers$Esri.WorldImagery,group="Esri",layerId="Esri") %>%
#       addRasterImage(I_ra,colors=col1,method="ngb",group = "Observed qual") %>%
#       addLegend(pal =  col1, values=IST$BIOIDX_TXG,group = "Observed qual", position = "topright",opacity=0.8,title="Observed quality")
# 
#     ### Finalize leaflet
#     m3 <- m %>%
#       addLayersControl(baseGroups= c("Positron","Esri"),overlayGroups = c("Observed qual"),position="topleft")
#     m3
#     mapview::mapshot(m3, url =paste0("D:/SIG/SIG_BAFU/leaflets/",GUILD,"_plausi.html"))
#   }
# 
# 
# # Leaflets plausibilisation  -------------------------------------------------------
# 
#   pacman::p_load(leafgl,leaflet,leafem,leafpop,data.table,raster,tidyverse,sf)
#   res <- LOAD("ha2add","data",dir)
#   st_geometry(res) <- "geometry"
# 
#   res2 <- st_transform(res,2056)
#   write_sf(obj = res2, dsn = paste0("D://SIG/SIG_BAFU/GPKG/G",GUILD,"_EPSG2056.gpkg"),layer="BSS")
#   rgdal::ogrListLayers(paste0("D://SIG/SIG_BAFU/GPKG/G",GUILD,".gpkg"))
# 
#   res2 <- st_transform(res,4326)
#   V1 <- col.bin(res2,bench,min.size=5,quiet=F)
# 
#   summary(V1[[1]]$potential_qual-V1[[1]]$Erganzungsbedarf)
# 
#   m = leaflet() %>%
#     addProviderTiles(provider = providers$CartoDB.Positron,group="Positron",layerId="Positron") %>%
#     addProviderTiles(provider = providers$Esri.WorldImagery,group="Esri",layerId="Esri") %>%
#     addPolygons(data = V1[[1]], group = "Erganzungsbedarf",weight=1,color="grey",opacity=0.8,fillColor = ~col,fillOpacity = 0.8,highlightOptions = highlightOptions(color = "black", weight = 2),popup = popupTable(V1[[1]][,c("BV.id","sp.in.BV","sp.in.bench","observed_qual","potential_qual","Erganzungsbedarf")],feature.id = FALSE,row.numbers = FALSE)) %>%
#     addLegend(colors=V1[[2]]$col,labels=V1[[2]]$lab,group = "Erganzungsbedarf", position = "topright",opacity=0.8,title=paste0("Erganzungsbedarf [ha] (max:",max(V1[[1]]$Erganzungsbedarf,na.rm=T),")")) %>%
#     addPolygons(data = V1[[1]], group = "BV",weight=1,color="grey",opacity=0.8,fillColor = "transparent",fillOpacity = 0,highlightOptions = highlightOptions(color = "black", weight = 2),popup = popupTable(V1[[1]][,c("BV.id","sp.in.BV","sp.in.bench","observed_qual","potential_qual","Erganzungsbedarf")]))
# 
# 
#   m$dependencies = c(m$dependencies,
#                      mapview:::popupLayoutDependencies())
# 
#   ## ADD IST
#   PRED <- LOAD("qual4leaf","data",dir)
#   st_geometry(PRED) <- "geometry"
#   IST <- st_transform(PRED[!is.na(PRED$ist),],3857)
#   grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
#   I_ra <- fasterize::fasterize(IST,grid2,field="ist")
#   col1 <- colorNumeric(palette = "viridis",IST$ist,na.color="#00000000")
#   # raster::writeRaster(I_ra,paste0("shp/",GUILD,"_IST.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
# 
#   m1 <- m %>%
#     addRasterImage(I_ra,colors=col1,method="ngb",group = "Observed qual") %>%
#     addLegend(pal =  col1, values=IST$ist,group = "Observed qual", position = "topright",opacity=0.8,title="Observed quality")
# 
# 
#   ## ADD SOLL
#   P <- setDT(PRED)[quality=="pred",]
#   st_geometry(P) <- "geometry"
#   P <- st_transform(P,3857)
#   grid2 <- raster::raster("D://SIG/data/grid100/grid100_3857.tif")
#   P_ra <- fasterize::fasterize(P,grid2,field="consensus")  ## replace by "consensus" when priorisation is done
#   # raster::writeRaster(P_ra,paste0("shp/",GUILD,"_SOLL2.tif"),format = 'GTiff', options=c("COMPRESS=DEFLATE", "PREDICTOR=2"),overwrite=TRUE)
# 
#   # col3 = colourvalues::colour_values_rgb(P$prio, palette = "inferno",include_alpha = FALSE)
#   col3 <- colorNumeric(palette = "inferno",P$consensus,na.color="#00000000",reverse=T)
# 
#   m2 <- m1 %>%
#     addRasterImage(P_ra,colors=col3,method="ngb",group = "Potential qual",maxBytes=4454211) %>%
#     addLegend(color =  col3(1),label="" ,group= "Potential qual", position = "topright",opacity=0.8,title="Potential quality")
# 
# 
#   ### Finalize leaflet
#   m3 <- m2 %>%
#     addLayersControl(baseGroups= c("Positron","Esri"),overlayGroups = c("Erganzungsbedarf","BV","Observed qual","Potential qual"),position="topleft") %>%
#     hideGroup(c("Observed qual","Potential qual"))
#   m3
#   mapview::mapshot(m3, url =paste0("D:/SIG/SIG_BAFU/leaflets/",GUILD,"_plausi.html"))
#   }
# 
# 

