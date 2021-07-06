pacman::p_load(sf, data.table,tidyverse)
source('~/R/infraeco_2021/code/functions_OI.R', encoding = 'UTF-8')

#' @title BSS_summary
#' @author Ervan Rutishauser
#' @description a wrapper to compute BSS by shape (provided by the user)
#' @param shape is a shapefile where a summary stat is returned for every shape 
#' @param shape.ID, shape identifier to be displayed in the summary table (e.g. Canton's name). By default the first column is used.
#' @return a table summarizing BSS by shape
#' @export
#' @examples
#' # default 
#' 
#' 

BSS_summary <- function(shape,shape.ID=NULL) {

  # Creation table BSS over all guilds --------------------------------------
  # LOCATION <- "C://Dossier_Ervan/R/"# location of all folders (IF office)
  # G <- c(2:8,10,12:17,19:20,22,25,26)
  # a=0
  # for (i in c(2:8,10,12:17,19:20,22,25,26)) {
  #   a=a+1
  #   dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",i,"-"),list.files(LOCATION))])
  #   guild <- readxl::read_xlsx("C:/Dossier_Ervan/R/infraeco_2021/data/Info_guildes.xlsx",sheet=1,col_names = T,col_types = "guess")
  #   bench <- guild[which(guild$ID==i),]$bench
  #   guild.name <- guild[which(guild$ID==i),]$name_guild_de
  #   
  #   res <- LOAD("ha2add","data",dir)
  #   res <- col.bin(res,bench,min.size=5)
  #   res <- res[[1]]
  #   
  #   ## Create table with all BSS
  #   if (i==2){
  #     TAB <- res[,c("BV.id","Erganzungsbedarf")]
  #     names(TAB) <- c("BV.id",i,"geometry")
  #   } else {
  #     TAB <- merge(TAB,res[,c("BV.id","Erganzungsbedarf")] %>% st_drop_geometry())
  #     names(TAB)[a+1] <- i
  #   }
  # }
  # save(TAB,file="C:/Dossier_Ervan/R/Miscallenious/data/BSS_all_guilds.Rdata")
  load("C:/Dossier_Ervan/R/Miscallenious/data/BSS_all_guilds.Rdata")  ##data = TAB BSS over all guilds

  
  ### Merge grid with shape
  # shape <- "D:/SIG/data/Limites/Limites_communales/swissBOUNDARIES3D_1_3_TLM_HOHEITSGEBIET.shp"
  # shape.id <- "NAME" 
  shape <- st_read(shape)
  
  # Convert into EPSG 2056
  load("C:/Dossier_Ervan/R/Grid100/grid100_sf_2056.Rdata")
  
  st_geometry(grid_sf) <- "centro"
  if(!identical(st_crs(grid_sf),st_crs(shape))) {
    shape <- st_transform(shape,st_crs(grid_sf))
  }
  grid <- st_join(grid_sf,shape[,c(shape.id)],join=st_intersects)
  

  ## Number of shape in BVs
  B <- grid_sf[,length(unique(shape.id)),by=BV04]

  ## Nbr of ha in each BV by guild
  load("C:/Dossier_Ervan/R/Miscallenious/data/data_final_stat.Rdata")  # data = G
  IDX <- match(c("CNHA",paste0("P",guild)),names(G))
  TT <- merge(setDT(G)[,..IDX],setDT(grid)[,c("CNHA","BV04",shape.id),with=FALSE],by="CNHA",all.x=TRUE)
  TT2 <- TT[!is.na(shape.id),]  # remove non-assigned pixel

  A <- TT2 %>% 
    dplyr::select(P2:BV04) %>% 
    group_by(BV04) %>% 
    summarise_all(list(function(x) sum(!is.na(x))))


  ## Number of ha in shapes by BV 
  idx <- match(c(shape.id,"BV04"),names(TT2))
  B <- TT2 %>% 
    group_by_at(idx) %>% 
    summarise_all(list(function(x) sum(!is.na(x))))


  ## proportion of quality in shape by BV
  D <- cbind(B[,1:2],B[, -c(1:3)]/A[match(B$BV04, A$BV04),-1])

  # Calcul des BSS par BV et shape
  col.id <- substr(names(D[, -c(1:2)]),2,3)
  IDX <- match(B$BV04, TAB$BV.id)
  E <- cbind(D[,c(1:2)],round(D[, -c(1:2)]*as_tibble(setDT(TAB)[IDX,..col.id])))

  ## Sum BSS by shape
  
  EE <- E %>% 
    group_by_at(1) %>% 
    summarise(across(names(D[, -c(1:2)]), sum,na.rm=TRUE))

  EE$total <- apply(EE[,2:ncol(EE)],1,sum,na.rm=TRUE)
  library(DT)
  datatable(EE)
  # fwrite(EE,file="C:/Dossier_Ervan/R/Miscallenious/data/BSS_by_communes.csv",sep=";")
} 


#' @title BSS_selection
#' @author Ervan Rutishauser
#' @description a wrapper to select BSS by shape (provided by the user). 
#' @param TAB , a summary table done through BSS_summary 
#' @param guild.id , a vector (numeric) with the desired guilds (c(2,4,7))
#' @param prio , indices of priorisation to be used ("guild_overlap","connectivity","historic_quality","env_suitability"). By default, all indices are used (and summed up).
#' @return a table summarizing BSS by shape
#' @export
#' @examples
#' # default 
#' 
#' 
BSS_selection <- function(shape,shape.id=NULL,guild.id=NULL,prio=NULL)
  # defin selected guilds
  if(!is.null(guild.id)){
    guild.id 
  } else {
    guild.id <- c(2:8,10,12:17,19,20,22,25,26)
  }

  # select index of priorisation
  load("C:/Dossier_Ervan/R/Grid100/grid100_sf_2056.Rdata")  # data = grid
  LOCATION <- "C://Dossier_Ervan/R/"# location of all folders (IF office)
  source("C:/Dossier_Ervan/R/infraeco_2021/code/functions_OI.R")
  
 
  
  if(is.null(prio)){
    prio <- "consensus" 
  } else {
    prio
  }
  
  GG <- setDT(grid)[,c("CNHA")]
  
  a=1
  for (i in guild.id) {
      a=a+1
      dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",i,"-"),list.files(LOCATION))])
      PRIO <- LOAD("prio","data",dir,silent=TRUE)
      GG <- merge(GG,setDT(PRIO)[,c("CNHA",prio),with=FALSE],by="CNHA",all.x=TRUE)
      names(GG)[a] <- paste0("P",i)
    }

  cols <- paste0("P",guild.id)
  
  over <- GG # create a new data set to keep info on quality index
  over[,cols] <- GG[,replace(.SD, !is.na(.SD), 1),.SDcols=cols]
  over[,"overlap" := ifelse(rowSums(.SD, na.rm=T)==0,NA,rowSums(.SD, na.rm=T)),.SDcols=cols]
  GG$overlap <- over$overlap
  
  ### Need to assign shape.id from shape to grid_sf here !!!!
  .
  .
  .
  .
  
  
  ## Nbr of ha Qp in each shape by guild
  IDX <- match(c("CNHA",paste0("P",guild.id)),names(GG))
  TT <- merge(setDT(GG)[,..IDX],setDT(grid)[!is.na(shape.id),c("CNHA","BV04",shape.id),with=FALSE],by="CNHA",all.x=TRUE)  
  TT2 <- TT[!is.na(shape.id),]  # remove non-assigned pixel
  
  idx <- match(c(shape.id),names(TT2))
  B <- TT2 %>% 
    # select(-any_of("CNHA"))  %>% 
    group_by_at(idx) %>% 
    summarise_all(list(function(x) sum(!is.na(x))))

  EE <- fread("C:/Dossier_Ervan/R/Miscallenious/data/BSS_by_canton.csv")
  names(EE)[1] <- shape.id
  GG <- merge(setDT(GG),setDT(grid)[!is.na(shape.id),c("CNHA",shape.id),with=FALSE],by="CNHA",all.x=TRUE)  
  over <- merge(setDT(over),setDT(grid)[!is.na(shape.id),c("CNHA",shape.id),with=FALSE],by="CNHA",all.x=TRUE)  
   
  # G[EE, on=.(canton), by=.EACHI,
   #     .(.(mapply(function(x, n) -head(sort(-x, partial=n), n),
   #                x=mget(cols), n=mget(paste0("i.", cols)), SIMPLIFY=FALSE)))]
  
  library(tidyverse)
  
  ## Select CNHA with lowest overlap (max hectares)
  get_CNHA_max <- function(.x, .group, .n_max, .dat) {
    .dat %>% 
      filter(canton == .group) %>% 
      filter(!is.na(.data[[.x]])) %>%
      # distinct(.data[[.x]],.keep_all=TRUE) %>% 
      arrange(.data$overlap,desc(.data[[.x]])) %>%
      slice(seq_len(.n_max)) %>% 
      pull(.data$CNHA)
  }
  
  Smax <- as_tibble(EE) %>%
    dplyr::select(-total) %>%
    pivot_longer(-c(canton)) %>%
    rowwise() %>% 
    mutate(CNHA = list(
      get_CNHA_max(name, canton, value, GG)
    )) %>%
    ungroup()
  
  # ## Check
  # GG[GG$CNHA%in%dat$CNHA[[1]],]
  # table(GG$overlap[GG$CNHA%in%Smin$CNHA[[1]]])
  # table(GG$overlap[!GG$CNHA%in%dat$CNHA[[1]] & GG$canton=="Aargau" & !is.na(GG$P2)])
  # GG[ GG$canton=="Aargau" & GG$overlap==1 & !is.na(GG$P2),]
  # table(duplicated(unlist(Smin$CNHA)))
  
  # remove duplicates
  un <- unlist(Smax$CNHA)
  res1 <- Map(`[`, Smax$CNHA, relist(!duplicated(un), skeleton = Smax$CNHA))
  # table(duplicated(unlist(res)))
  # table(GG$overlap[GG$CNHA%in%res[[1]]])
  A <- unlist(lapply(Smin$CNHA,length)) - unlist(lapply(res,length))
  head(A)
  

  ## Select CNHA with maximum overlap (min area)
  get_CNHA_min <- function(.x, .group, .n_max, .dat) {
    .dat %>% 
      filter(canton == .group) %>% 
      filter(!is.na(.data[[.x]])) %>%
      # distinct(.data[[.x]],.keep_all=TRUE) %>% 
      arrange(desc(.data$overlap),desc(.data[[.x]])) %>%   ## filter overlap by descending order
      slice(seq_len(.n_max)) %>% 
      pull(.data$CNHA)
  }
  

  Smin <- as_tibble(EE) %>%
    dplyr::select(-total) %>%
    pivot_longer(-c(canton)) %>%
    rowwise() %>% 
    mutate(CNHA = list(
      get_CNHA_min(name, canton, value, GG)
    )) %>%
    ungroup()
  
  # remove duplicates
  un <- unlist(Smin$CNHA)
  res2 <- Map(`[`, Smin$CNHA, relist(!duplicated(un), skeleton = Smin$CNHA))
 
  
  # Map results
  load("C:/Dossier_Ervan/R/Grid100/grid100_sf.Rdata")
  G_sf <- merge(setDT(GG),setDT(grid_sf)[,c("CNHA","geometry"),with=FALSE],by="CNHA",all.x=TRUE)  
  st_geometry( G_sf ) <- "geometry"
 
  MIN <- G_sf[G_sf$CNHA%in%res1[[1]],]
  summary(MIN$P2)
  MIN <- st_transform(MIN,4326)
  MAX <- G_sf[G_sf$CNHA%in%res2[[1]],]
  summary(MAX$P2)
  MAX <- st_transform(MAX,4326)
  
  pacman::p_load(leafgl,leaflet,leafem,leafpop,data.table,raster,sf,tidyverse)
  pal <- colorNumeric("viridis", MIN$P2)
  pal2 <- colorNumeric("plasma", MAX$P2)
  
  m <- leaflet() %>%
        addProviderTiles(provider = providers$CartoDB.Positron,group="Positron",layerId="Positron") %>%
        addProviderTiles(provider = providers$Esri.WorldImagery,group="Esri",layerId="Esri") %>%
       addPolygons(data = MIN, group = "min",weight=1,color="grey",opacity=0.8,fillColor = pal(MIN$P2),fillOpacity = 0.8,highlightOptions = highlightOptions(color = "black", weight = 2),popup = popupTable(MIN[,c("CNHA","P2","overlap")],feature.id = FALSE,row.numbers = FALSE)) %>%
      addLegend(pal = pal,values=MIN$P2,group = "min", position = "topright",opacity=0.8,title="BSS min") %>%
      addPolygons(data = MAX, group = "max",weight=1,color="grey",opacity=0.8,fillColor = pal2(MAX$P2),fillOpacity = 0.8,highlightOptions = highlightOptions(color = "black", weight = 2),popup = popupTable(MAX[,c("CNHA","P2","overlap")],feature.id = FALSE,row.numbers = FALSE)) %>%
    addLegend(pal = pal2, values = MAX$P2,group = "max", position = "topright",opacity=0.8,title="BSS MAX")
 
    
  m <-  m %>%
        addLayersControl(baseGroups= c("Positron","Esri"),overlayGroups = c("min","max"),position="topleft")
      
  m
 
  