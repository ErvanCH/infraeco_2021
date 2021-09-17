#### Summary statistics among all guilds
library(raster)
library(sf)
library(data.table)

# Define protected area ---------------------------------------------------------
# #### Compile information on protected area
# load("C:/Dossier_Ervan/R/Grid100/grid100_sf_with_enviro.Rdata")
# st_geometry(grid_sf) <- "centro"

# ## Zone de protection 1: Inventaires/réserves forestières et PN /Sites Ramsar et Emeraude
# TARG <- c("D://SIG/data/Aires_protegees/Reserve_ProNat/NSG_PN_2017_lv95.shp","D:/SIG/data/Aires_protegees/Reserve_ProNat/Wald_PN_2017_lv95.shp","D:/SIG/data/Aires_protegees/Emeraude_LV95/smaragd.shp","D://SIG/data/Aires_protegees/Waldreserven/waldreservate.shp","D:/SIG/data/Parcs_suisse_LV95/N2020_Revision_Park_20200101.shp")
# 
# 
# for (i in 1:length(TARG)) {
# r1 <- st_read(TARG[i])
# if (i==5) {
#   r1 <- r1[r1$Kategorie=="SNP",]  ## ne garde que le parc national
# }
# r1$VV <- 1
# r1 <- st_transform(r1,st_crs(grid_sf))
# system.time(grid_sf <- st_join(grid_sf,r1[,"VV"],join=st_intersects))
# names(grid_sf)[grep("VV",names(grid_sf))] <- paste0("AP_1_",i)
# }
# 
# 
# ## Zone de protection 2: district francs et Parcs Naturels
# TARG2 <- c("D://SIG/data/Aires_protegees/Jagdbanngebiete_LV95/JB.shp","D:/SIG/data/Parcs_suisse_LV95/N2020_Revision_Park_20200101.shp")
# 
# for (i in 1:length(TARG2)) {
#   r1 <- st_read(TARG2[i])
#   if (i==2) {
#     r1 <- r1[r1$Kategorie!="SNP",]  ## enlève le parc national
#   }
#   r1$VV <- 1
#   r1 <- st_transform(r1,st_crs(grid_sf))
#   system.time(grid_sf <- st_join(grid_sf,r1[,"VV"],join=st_intersects))
#   names(grid_sf)[grep("VV",names(grid_sf))] <- paste0("AP_2_",i)
# }
# 
# G <- grid_sf[,c("grid.id", "BV04", "CX", "CY", "CNHA", "bioregion", "forest", 
#                 "urban", "lac", "subreg", "amphi", "au", "flachm", "hochm", 
#                 "TWW", "vogel","AP_1_1","AP_1_2", "AP_1_3", "AP_1_4", "AP_2_1", "AP_2_2", "AP_1_5", "geometry", 
#                 "centro")]
# G[,"prot1" := ifelse(au==1|flachm==1|hochm==1|TWW==1|vogel==1|amphi==1|AP_1_1|AP_1_2|AP_1_3|AP_1_4|AP_1_5,1,0)]
# G[,"prot2" := ifelse(AP_2_1|AP_2_2,1,0)]
# 
# save(G,file="C:/Dossier_Ervan/R/Grid100/grid100_sf_with_protection.Rdata")

# Qualité existante -------------------------------------------------------

load("C:/Dossier_Ervan/R/Grid100/grid100_sf_with_protection.Rdata")
GUILD <- c(2:10,12:17,19:20,22,24,25,26)  # manque 1,11,18,21 / IST only: 9 et 24


a=ncol(G)

for (i in 1:length(GUILD)) {
  a=a+1
  FF <- grep(paste0("_",GUILD[i],"_bio-idx",collapse="|"),list.files("C:/Dossier_Ervan/R/infraeco_2021/data/infofauna"),value=T)
  ist <- fread(paste("C:/Dossier_Ervan/R/infraeco_2021/data/infofauna", FF,sep="/"))
  if (GUILD[i] %in% c(24,25,26)) {
    G <- merge(setDT(G),setDT(ist)[,c("CNHA","tx_ALL_40")],by="CNHA",all.x=T)
    names(G)[a] <- paste0("G",GUILD[i])
  } else {
    G <- merge(setDT(G),setDT(ist)[,c("CNHA","BIOIDX_TXG")],by="CNHA",all.x=T)
    names(G)[a] <- paste0("G",unique(ist$GUILD))
  }
}

G <- setDT(G)[!duplicated(CNHA),]
cols <- grep("G",names(G))

QUAL <- G[, lapply(.SD, function(x) ifelse(x< median(x,na.rm=T),1,2)), .SDcols = cols]  # binarize quality
QUAL[,"Qhigh" := apply(.SD,1, function(x) any(x==2))]
QUAL[is.na(Qhigh),"Qlow" := apply(.SD,1, function(x) any(x==1))]

G$Qhigh <- QUAL$Qhigh
G$Qlow <- QUAL$Qlow

G[,cols] <- G[,replace(.SD, !is.na(.SD), 1),.SDcols=cols]
G[,"overlap" := ifelse(rowSums(.SD, na.rm=T)==0,NA,rowSums(.SD, na.rm=T)),.SDcols=cols]
hist(G$overlap)

## Nbr. d'hectare unique
B <- setDT(G)[,.("n.obs"=length(which(!is.na(overlap))),"ntot"=.N)]
B[,"prop":=round((n.obs*100)/ntot,2)]
B

table(G[!is.na(overlap),prot1]) # Qualité en protection 1
table(G[!is.na(overlap) & is.na(prot1),prot2]) # Qualité en protection 2
G[!is.na(overlap) & is.na(prot1) & is.na(prot2),.N] # Qualité sans protection

nrow(G[!is.na(Qhigh)])
table(G[!is.na(Qhigh),prot1]) # Qualité en protection 1
table(G[!is.na(Qhigh) & is.na(prot1),prot2]) # Qualité en protection 2
G[!is.na(Qhigh) & is.na(prot1) & is.na(prot2),.N] # Qualité sans protection

nrow(G[!is.na(Qlow)])
table(G[!is.na(Qlow),prot1]) # Qualité en protection 1
table(G[!is.na(Qlow) & is.na(prot1),prot2]) # Qualité en protection 2
G[!is.na(Qlow) & is.na(prot1) & is.na(prot2),.N] # Qualité sans protection

protected<- c(length(G$prot1),length(na.exclude(G$prot1)),length(na.exclude(G$prot2)),length(G$prot1)-(length(na.exclude(G$prot1))+length(na.exclude(G$prot2))))
total.quality.protected<-c(length(na.exclude(G$Qhigh))+length(na.exclude(G$Qlow)),table(G[!is.na(Qhigh) | !is.na(Qlow),prot1]), table(G[!is.na(Qhigh) | !is.na(Qlow),prot2]),length(na.exclude(G$Qhigh))+length(na.exclude(G$Qlow))-(table(G[!is.na(Qhigh) | !is.na(Qlow),prot1]) + table(G[!is.na(Qhigh) | !is.na(Qlow),prot2])))
quality1.protected<-c(length(na.exclude(G$Qlow)),table(G[!is.na(Qlow),prot1]), table(G[!is.na(Qlow),prot2]),length(na.exclude(G$Qlow))-(table(G[!is.na(Qlow),prot1]) + table(G[!is.na(Qlow),prot2])))
quality2.protected<-c(length(na.exclude(G$Qhigh)),table(G[!is.na(Qhigh),prot1]), table(G[!is.na(Qhigh),prot2]),length(na.exclude(G$Qhigh))-(table(G[!is.na(Qhigh),prot1]) + table(G[!is.na(Qhigh),prot2])))
no.quality.protected<-protected-total.quality.protected
quality.protection<-rbind(protected,total.quality.protected,quality1.protected,quality2.protected, no.quality.protected )
colnames(quality.protection)<-c('total', 'strict_protection','light_protection','no_protection')
row.names(quality.protection)<-c('ha_in_CH', 'ha_with_quality','ha_quality_level1','ha_quality_level2','ha_no_quality')
write.csv2(quality.protection, file = 'C:/Dossier_Ervan/R/Miscallenious/data/quality_in_protection.csv',quote =F)

# Hectares to add ---------------------------------------------------------
LOCATION <- "C://Dossier_Ervan/R/"# location of all folders (IF office)
G <- c(2:8,10,12:17,19:20,22,25,26)
TAB <- matrix(NA,length(G),2)
a=0
for (i in c(2:8,10,12:17,19:20,22,25,26)) {
  a=a+1
  dir <- paste0(LOCATION,list.files(LOCATION)[grep(paste0("G",i,"-"),list.files(LOCATION))])
  guild <- readxl::read_xlsx("data/Info_guildes.xlsx",sheet=1,col_names = T,col_types = "guess")
  bench <- guild[which(guild$ID==i),]$bench
  guild.name <- guild[which(guild$ID==i),]$name_guild_de
  
  res <- LOAD("ha2add","data",dir)
  res <- col.bin(res,bench,min.size=5)
  res <- res[[1]]
  
  ## Create table with all BSS
  if (i==2){
    temp <- res[,c("BV.id","Erganzungsbedarf")]
    names(temp) <- c("BV.id",i,"geometry")
  } else {
    temp <- merge(temp,res[,c("BV.id","Erganzungsbedarf")] %>% st_drop_geometry())
    names(temp)[a+1] <- i
  }
  
  # sum BSS
  TAB[a,2] <- sum(res$Erganzungsbedarf,na.rm=T)
  
}
TAB[,1] <- c(2:8,10,12:17,19:20,22,25,26)
sum(TAB[,2])/B$ntot
