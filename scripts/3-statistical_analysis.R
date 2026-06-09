setwd("~/GitHub/PhIr_Pfrac")
source("scripts/0-packages.R")

soils_PSI <- read.csv("formatted spreadsheets/soils_PSI.csv")
thaw <- read.csv("formatted spreadsheets/thaw_sample.csv")
horizon_depth <- read.csv("raw data/core_lengths.csv")
horizon_depth <- horizon_depth %>% dplyr::select(SampleEvent,Area,SiteName,PlotNumber,Org_end) %>% rename(depth = Org_end)

soils_PSI <- merge(soils_PSI, horizon_depth, by=intersect(names(soils_PSI[c(1,2,3,5)]), names(horizon_depth[,1:4])))
soils_PSI <- soils_PSI %>% mutate(SiteName= factor(SiteName, levels = c("Dry", "Mesic", "Hydric"))) %>% mutate(Area=factor(Area)) %>% mutate(SampleEvent=factor(SampleEvent,levels = c("Early Season", "Mid Season", "Late Season"))) %>% mutate(treatment=factor(treatment,levels=c("Acidic Dry", "Acidic Mesic", "Acidic Hydric", "Non-acidic Dry", "Non-acidic Mesic", "Non-acidic Hydric")))

thaw<- thaw %>% mutate(Site= factor(Site, levels = c("Dry", "Mesic", "Hydric")))%>% mutate(Area=factor(Area)) %>% mutate(SampleEvent=factor(SampleEvent, levels = c("Early Season", "Mid Season", "Late Season"))) %>% mutate(treatment = paste(Area,Site,sep = " "))

thaw_summ <- thaw %>% dplyr::select(SampleEvent,treatment,Value) %>% group_by(SampleEvent,treatment) %>%
  summarize(mean = mean(Value, na.rm = TRUE),
            sd=sd(Value, na.rm = TRUE),
            se=sd(Value)/sqrt(n()),
            n =n())



soils_PSI$lnH2OSRP <- log(soils_PSI$PfracH2OSRP)
shapiro.test(soils_PSI$lnH2OSRP)
soils_PSI$lnHClP <- log(soils_PSI$PfracHClP)
shapiro.test(soils_PSI$lnHClP)
soils_PSI$lnPSI <- log(soils_PSI$PfracPSI)
shapiro.test(soils_PSI$lnPSI)
soils_PSI$H2Oinv <- 1/soils_PSI$PfracH2OSRP
soils_PSI$sq_Fe <- sqrt(soils_PSI$NonCFe)
shapiro.test(soils_PSI$sq_Fe)

### summary tables ####
table_stacked <- soils_PSI %>% dplyr::select(treatment,SampleEvent,moisture,depth,PfracPSI,pH,LOI,C.,N.,CN_ratio,PfracBDFe,PfracBDSRP,PfracH2OSRP,PfracHClCa,PfracHClP,PfracNaOHAl,PfracNaOHSRP,PfracNRP,PfracRP,PfracTotAl,PfracTotCa,PfracTotFe,PfracTotP,ResP,FefracDHFe,OrgFe,NonCFe,ResFe) %>% gather("Fraction", "Value", 3:28)

table_summary_time <- table_stacked %>% group_by(SampleEvent,treatment,Fraction) %>% 
  summarize(mean=mean(Value, na.rm = TRUE),
            sd=sd(Value, na.rm = TRUE),
            se = sd(Value)/sqrt(n()),
            n=length(Value),)
table_summary <- table_stacked %>% group_by(treatment,Fraction) %>% 
  summarize(mean=mean(Value, na.rm = TRUE),
            sd=sd(Value, na.rm = TRUE),
            se=sd(Value)/sqrt(n()),
            n=length(Value))


summary_time <-function(raw,value,sheet){
  if (sheet==0) {
    AV<- raw[which(raw$Fraction == value),c(1,2,4)] %>% spread(treatment,mean)
    SE<- raw[which(raw$Fraction == value),c(1,2,6)] %>% spread(treatment,se)
    AV[,2:7] <- round(AV[,2:7],1)
    SE[,2:7] <- round(SE[,2:7],1)
  } else if (sheet==1) {
    AV<- raw[,c(1,2,3)] %>% spread(treatment,mean)
    SE<- raw[,c(1,2,5)] %>% spread(treatment,se)
    AV[,2:7] <- round(AV[,2:7],1)
    SE[,2:7] <- round(SE[,2:7],1)
  }
  stats_summary <- data.frame(SampleEvent = c("Early Season", "Mid Season","Late Season"))
  stats_summary$'Acidic Dry' <- paste(AV$'Acidic Dry',"±",SE$'Acidic Dry')
  stats_summary$'Acidic Mesic' <- paste(AV$'Acidic Mesic',"±",SE$'Acidic Mesic')
  stats_summary$'Acidic Hydric' <- paste(AV$'Acidic Hydric',"±",SE$'Acidic Hydric')
  stats_summary$'Non-acidic Dry' <- paste(AV$'Non-acidic Dry',"±",SE$'Non-acidic Dry')
  stats_summary$'Non-acidic Mesic' <- paste(AV$'Non-acidic Mesic',"±",SE$'Non-acidic Mesic')
  stats_summary$'Non-acidic Hydric' <- paste(AV$'Non-acidic Hydric',"±",SE$'Non-acidic Hydric')
  stats_summary$SampleEvent <- factor(stats_summary$SampleEvent, levels = c("Early Season","Mid Season","Late Season"))
  stats_summary<- stats_summary[order(stats_summary$SampleEvent),]
  return(stats_summary)
}

summary <-function(raw){
  AV<- raw[,c(1,2,3)] %>% spread(Fraction, mean)
  SE<- raw[,c(1,2,5)] %>% spread(Fraction, se)
  AV[,2:26] <- round(AV[,2:26],1)
  SE[,2:26] <- round(SE[,2:26],1)
  stats_summary <- data.frame(treatment = c("Acidic Dry", "Acidic Mesic","Acidic Hydric","Non-acidic Dry","Non-acidic Mesic","Non-acidic Hydric"))
  stats_summary$'depths' <- paste(AV$'depth',"±",SE$'depth')
  stats_summary$'LOI' <- paste(AV$'LOI',"±",SE$'LOI')
  stats_summary$'pH' <- paste(AV$'pH',"±",SE$'pH')
  stats_summary$'C.' <- paste(AV$'C.',"±",SE$'C.')
  stats_summary$'N.' <- paste(AV$'N.',"±",SE$'N.')
  stats_summary$'CN_ratio' <- paste(AV$'CN_ratio',"±",SE$'CN_ratio')
  stats_summary$'PfracNRP' <- paste(AV$'PfracNRP',"±",SE$'PfracNRP')
  stats_summary$'PfracH2OSRP' <- paste(AV$'PfracH2OSRP',"±",SE$'PfracH2OSRP')
  stats_summary$'PfracBDSRP' <- paste(AV$'PfracBDSRP',"±",SE$'PfracBDSRP')
  stats_summary$'PfracNaOHSRP' <- paste(AV$'PfracNaOHSRP',"±",SE$'PfracNaOHSRP')
  stats_summary$'PfracHClP' <- paste(AV$'PfracHClP',"±",SE$'PfracHClP')
  stats_summary$'PfracRP' <- paste(AV$'PfracRP',"±",SE$'PfracRP')
  stats_summary$'ResP' <- paste(AV$'ResP',"±",SE$'ResP')
  stats_summary$'PfracTotP' <- paste(AV$'PfracTotP',"±",SE$'PfracTotP')
  stats_summary$'PfracBDFe' <- paste(AV$'PfracBDFe',"±",SE$'PfracBDFe')
  stats_summary$'OrgFe' <- paste(AV$'OrgFe',"±",SE$'OrgFe')
  stats_summary$'NonCFe' <- paste(AV$'NonCFe',"±",SE$'NonCFe')
  stats_summary$'FefracDHFe' <- paste(AV$'FefracDHFe',"±",SE$'FefracDHFe')
  stats_summary$'ResFe' <- paste(AV$'ResFe',"±",SE$'ResFe')
  stats_summary$'PfracTotFe' <- paste(AV$'PfracTotFe',"±",SE$'PfracTotFe')
  stats_summary$'PfracNaOHAl' <- paste(AV$'PfracNaOHAl',"±",SE$'PfracNaOHAl')
  stats_summary$'PfracTotAl' <- paste(AV$'PfracTotAl',"±",SE$'PfracTotAl')
  stats_summary$'PfracHClCa' <- paste(AV$'PfracHClCa',"±",SE$'PfracHClCa')
  stats_summary$'PfracTotCa' <- paste(AV$'PfracTotCa',"±",SE$'PfracTotCa')
  stats_summary$'PfracPSI' <- paste(AV$'PfracPSI',"±",SE$'PfracPSI')
  stats_summary$treatment <- factor(stats_summary$treatment, levels = c("Acidic Dry", "Acidic Mesic","Acidic Hydric","Non-acidic Dry","Non-acidic Mesic","Non-acidic Hydric"))
  stats_summary<- stats_summary[order(stats_summary$treatment),]
  return(stats_summary)
}
moist_summ <- summary_time(raw = table_summary_time,value = "moisture", sheet = 0)
write.csv(moist_summ, "formatted spreadsheets/moisture_table.csv", row.names = F)
thaw_table <- summary_time(raw = thaw_summ,sheet = 1)
write.csv(thaw_table, "formatted spreadsheets/thaw_table.csv",row.names = F)

treatment_summ <- summary(raw = table_summary)

soil_conditions_table <- treatment_summ %>% dplyr::select(treatment,depths,pH,LOI,C.,N.,CN_ratio) %>% rename("Organic horizon depth (cm)"=depths,"Soil pH"=pH,"% organic matter"=LOI,"Carbon %"=C.,"Nitrogen %"=N.)
write.csv(soil_conditions_table,"formatted spreadsheets/soil_conditions_table.csv",row.names = F, fileEncoding = "UTF-8")

P_table <- treatment_summ %>% dplyr::select(treatment, PfracRP,PfracNRP, ResP,PfracTotP) %>% rename("Total rP"=PfracRP,"Total nrP"=PfracNRP,"Res~P"=ResP, "Total P"=PfracTotP)
write.csv(P_table,"formatted spreadsheets/TotalP_table.csv",row.names = F)

rP_table <- treatment_summ %>% dplyr::select(treatment,PfracH2OSRP,PfracBDSRP,PfracNaOHSRP,PfracHClP,PfracRP) %>% rename("loosely sorbed rP"=PfracH2OSRP,"rP~iron oxides"=PfracBDSRP,"rP~aluminum oxides"=PfracNaOHSRP,"P~calcareous minerals"=PfracHClP,"Total rP"=PfracRP)
write.csv(rP_table, "formatted spreadsheets/rP_table.csv",row.names = F)

iron_table <- treatment_summ %>% dplyr::select(treatment, PfracBDFe,OrgFe,NonCFe,FefracDHFe,ResFe,PfracTotFe) %>% rename("Iron oxides"=PfracBDFe,"Organic-bound iron" =OrgFe,"Non-crystalline iron"=NonCFe,"Crystalline iron"=FefracDHFe,"Residual iron"=ResFe, "Total Iron"=PfracTotFe)
write.csv(iron_table,"formatted spreadsheets/Iron_table.csv",row.names = F)

minerals_table <- treatment_summ %>% dplyr::select(treatment, PfracNaOHAl,PfracTotAl,PfracHClCa,PfracTotCa,PfracPSI) %>% rename("Aluminum oxides"=PfracNaOHAl, "Total Al"=PfracTotAl,"Calcium carbonate"=PfracHClCa,"Total Ca"=PfracTotCa, "Phosphate sorption index"=PfracPSI)
write.csv(minerals_table,"formatted spreadsheets/minerals_table.csv",row.names = F)

pvalue_table <- data.frame(parameters = c("TotalP","rP_Perc","nrP_Pec","resP_Perc","H2OrP","BDrP","NaOHrP","HClP","TotalFe","OrgFe","NonCFe","CFe","BDFe","NaOHAl","HClCa","TotalAl","TotalCa","PSI"), "soiltype"=NA,"hillslopeposition"=NA,"Hillslopeposition_Acidic"=NA,"Hillslopeposition_Non-acidic"=NA,"Growingseason"=NA)

pvalue_table_pt2 <- data.frame(parameters = c("H2OrP~BDrP","BDFe~ExtFe","H2OrP~OrgFe","H2OrP~BDFe","H2OrP~NCFe","H2OrP~CFe","H2OrP~NaOHAl_a","H2OrP~NaOHAl_na","HClP~HClCa_na","PSI~H2OrP"),"All_soils"=NA,"Hillslope"=NA)

pvalue_table_pt3 <- data.frame(parameters = c("H2OrP~BDrP","H2OrP~NaOHrP","H2OrP~HClP","BDrP~NCFe","BDrP~OrgFe","BDrP~CFe_a","H2OrP~HClCa_na","NaOHrP~NaOHAl_a","HClP~HClCa_na"), "All_soils" = NA,"Upland"=NA,"Midslope"=NA,"Lowland"=NA)


## Soil conditions ####
# soil pH
range(soils_PSI$pH[which(soils_PSI$Area=="Acidic")])
range(soils_PSI$pH[which(soils_PSI$Area=="Non-acidic")])

Area_pH <-lm(pH~Area, data = soils_PSI)
base::summary(Area_pH)
anova(Area_pH)
hist(resid(Area_pH))

# organic matter content
mean(soils_PSI$LOI[which(soils_PSI$SiteName=="Mesic")])
sd(soils_PSI$LOI[which(soils_PSI$SiteName=="Mesic")])/sqrt(length(soils_PSI$LOI[which(soils_PSI$SiteName=="Mesic")]))
mean(soils_PSI$LOI[which(soils_PSI$SiteName=="Dry")])
sd(soils_PSI$LOI[which(soils_PSI$SiteName=="Dry")])/sqrt(length(soils_PSI$LOI[which(soils_PSI$SiteName=="Dry")]))

mean(soils_PSI$LOI[which(soils_PSI$Area=="Acidic")])
sd(soils_PSI$LOI[which(soils_PSI$Area=="Acidic")])/sqrt(length(soils_PSI$LOI[which(soils_PSI$Area=="Acidic")]))
mean(soils_PSI$LOI[which(soils_PSI$Area=="Non-acidic")])
sd(soils_PSI$LOI[which(soils_PSI$Area=="Non-acidic")])/sqrt(length(soils_PSI$LOI[which(soils_PSI$Area=="Non-acidic")]))

LOI_site <- lm(LOI~SiteName, data = soils_PSI)
shapiro.test(resid(LOI_site))
anova(LOI_site)
base::summary(LOI_site)

# CN ratios
mean(soils_PSI$CN_ratio[soils_PSI$treatment=="Acidic Mesic"])
sd(soils_PSI$CN_ratio[soils_PSI$treatment=="Acidic Mesic"])/sqrt(length(soils_PSI$CN_ratio[soils_PSI$treatment=="Acidic Mesic"]))

mean(soils_PSI$CN_ratio[soils_PSI$treatment=="Non-acidic Hydric"])
sd(soils_PSI$CN_ratio[soils_PSI$treatment=="Non-acidic Hydric"])/sqrt(length(soils_PSI$CN_ratio[soils_PSI$treatment=="Non-acidic Hydric"]))

# moisture content
soils_PSI[which(soils_PSI$moisture==max(soils_PSI$moisture[which(soils_PSI$Area=="Acidic")])),c(1:4,10)]
soils_PSI[which(soils_PSI$moisture==max(soils_PSI$moisture[which(soils_PSI$Area=="Non-acidic")])),c(1:4,10)]
mean(soils_PSI$moisture[soils_PSI$treatment=="Acidic Hydric"&soils_PSI$SampleEvent=="Early Season"])
sd(soils_PSI$moisture[soils_PSI$treatment=="Acidic Hydric"&soils_PSI$SampleEvent=="Early Season"])/sqrt(length(soils_PSI$moisture[soils_PSI$treatment=="Acidic Hydric"&soils_PSI$SampleEvent=="Early Season"]))
mean(soils_PSI$moisture[soils_PSI$treatment=="Non-acidic Hydric"&soils_PSI$SampleEvent=="Early Season"])
sd(soils_PSI$moisture[soils_PSI$treatment=="Non-acidic Hydric"&soils_PSI$SampleEvent=="Early Season"])/sqrt(length(soils_PSI$moisture[soils_PSI$treatment=="Non-acidic Hydric"&soils_PSI$SampleEvent=="Early Season"]))


soils_PSI[which(soils_PSI$moisture==min(soils_PSI$moisture[which(soils_PSI$Area=="Acidic")])),c(1:4,10)]
soils_PSI[which(soils_PSI$moisture==min(soils_PSI$moisture[which(soils_PSI$Area=="Non-acidic")])),c(1:4,10)]
mean(soils_PSI$moisture[soils_PSI$treatment=="Acidic Hydric"&soils_PSI$SampleEvent=="Late Season"])
sd(soils_PSI$moisture[soils_PSI$treatment=="Acidic Hydric"&soils_PSI$SampleEvent=="Late Season"])/sqrt(length(soils_PSI$moisture[soils_PSI$treatment=="Acidic Hydric"&soils_PSI$SampleEvent=="Late Season"]))
mean(soils_PSI$moisture[soils_PSI$treatment=="Non-acidic Dry"&soils_PSI$SampleEvent=="Late Season"])
sd(soils_PSI$moisture[soils_PSI$treatment=="Non-acidic Dry"&soils_PSI$SampleEvent=="Late Season"])/sqrt(length(soils_PSI$moisture[soils_PSI$treatment=="Non-acidic Dry"&soils_PSI$SampleEvent=="Late Season"]))


### Total soil phosphorus ####
base::summary(soils_PSI$PfracTotP, na.rm = T)
TP_area <-lm(log(PfracTotP)~Area, data = soils_PSI)
pvalue_table[1,2] <-p_value_ext(TP_area, pt = 1)

TP_site <- lm(PfracTotP~SiteName, data = soils_PSI)
pvalue_table[1,3] <- p_value_ext(TP_site, pt = 1)

TP_time <- lm(PfracTotP~SampleEvent, data = soils_PSI)
pvalue_table[1,6] <- p_value_ext(TP_time, pt = 1)

# Total P in acidic soils
TP_a_site <- lm(PfracTotP~SiteName, data = soils_PSI[which(soils_PSI$Area == "Acidic"),])
pvalue_table[1,4] <- p_value_ext(TP_a_site, pt= 1)
  
mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Mesic")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Mesic")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Mesic")]))
mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Dry")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Dry")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Dry")]))
mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Hydric")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Hydric")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Hydric")]))

# Total P in non-acidic soils
TP_na_site <- lm(PfracTotP~SiteName, data = soils_PSI[which(soils_PSI$Area == "Non-acidic"),])
base::summary(TP_na_site)
pvalue_table[1,5] <- p_value_ext(TP_na_site, pt=1)

mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Hydric")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Hydric")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Hydric")]))
mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Dry")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Dry")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Dry")]))
mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Mesic")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Mesic")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Mesic")]))

# rP and total P
mean(soils_PSI$PfracRP/soils_PSI$PfracTotP, na.rm = TRUE)*100
sd(soils_PSI$PfracRP/soils_PSI$PfracTotP, na.rm = TRUE)/sqrt(length(soils_PSI$PfracRP))*100

rP_check <- soils_PSI %>% dplyr::select(SampleEvent, Area, SiteName, PlotNumber,PfracRP, PfracTotP)
rP_check$Perc <- rP_check$PfracRP/rP_check$PfracTotP*100
rP_perc_area <- lm(Perc~Area, data = rP_check)
pvalue_table[2,2] <- p_value_ext(rP_perc_area,pt=1)

rP_perc_site <- lm(Perc~SiteName, data = rP_check)
pvalue_table[2,3] <- p_value_ext(rP_perc_site,pt=1)

mean(soils_PSI$PfracRP[soils_PSI$SiteName=="Dry"]/soils_PSI$PfracTotP[soils_PSI$SiteName=="Dry"], na.rm = TRUE)*100
sd(soils_PSI$PfracRP[soils_PSI$SiteName=="Dry"]/soils_PSI$PfracTotP[soils_PSI$SiteName=="Dry"], na.rm = TRUE)/sqrt(length(soils_PSI$PfracRP[soils_PSI$SiteName=="Dry"]))*100

mean(soils_PSI$PfracRP[soils_PSI$SiteName=="Mesic"]/soils_PSI$PfracTotP[soils_PSI$SiteName=="Mesic"], na.rm = TRUE)*100
sd(soils_PSI$PfracRP[soils_PSI$SiteName=="Mesic"]/soils_PSI$PfracTotP[soils_PSI$SiteName=="Mesic"], na.rm = TRUE)/sqrt(length(soils_PSI$PfracRP[soils_PSI$SiteName=="Mesic"]))*100

mean(soils_PSI$PfracRP[soils_PSI$SiteName=="Hydric"]/soils_PSI$PfracTotP[soils_PSI$SiteName=="Hydric"], na.rm = TRUE)*100
sd(soils_PSI$PfracRP[soils_PSI$SiteName=="Hydric"]/soils_PSI$PfracTotP[soils_PSI$SiteName=="Hydric"], na.rm = TRUE)/sqrt(length(soils_PSI$PfracRP[soils_PSI$SiteName=="Hydric"]))*100

rP_perc_site_a <- lm(Perc~SiteName, data = rP_check[rP_check$Area=="Acidic",])
pvalue_table[2,4] <- p_value_ext(rP_perc_site_a,pt=1)

rP_perc_site_na <- lm(Perc~SiteName, data = rP_check[rP_check$Area=="Non-acidic",])
pvalue_table[2,5] <- p_value_ext(rP_perc_site_na,pt=1)

rP_perc_time <- lm(Perc~SampleEvent, data = rP_check)
pvalue_table[2,6] <- p_value_ext(rP_perc_time,pt=1)

# nrP and total P
mean(soils_PSI$PfracNRP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Acidic")], na.rm = TRUE)*100
sd(soils_PSI$PfracNRP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Acidic")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracNRP[which(soils_PSI$Area=="Acidic")]))*100
mean(soils_PSI$PfracNRP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Non-acidic")], na.rm = TRUE)*100
sd(soils_PSI$PfracNRP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Non-acidic")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracNRP[which(soils_PSI$Area=="Non-acidic")]))*100

NRP_check <- soils_PSI %>% dplyr::select(SampleEvent, Area, SiteName, PlotNumber,PfracNRP, PfracTotP)
NRP_check$Perc <- NRP_check$PfracNRP/NRP_check$PfracTotP*100
NRP_perc_area <- lm(Perc~Area, data = NRP_check)
pvalue_table[3,2] <- p_value_ext(NRP_perc_area,pt=1)

NRP_perc_site <- lm(Perc~SiteName, data = NRP_check)
pvalue_table[3,3] <- p_value_ext(NRP_perc_site, pt = 1)

NRP_perc_site_a <- lm(Perc~SiteName, data = NRP_check[NRP_check$Area=="Acidic",])
pvalue_table[3,4] <- p_value_ext(NRP_perc_site_a, pt =1)

NRP_perc_site_na <- lm(Perc~SiteName, data = NRP_check[NRP_check$Area=="Non-acidic",])
pvalue_table[3,5] <- p_value_ext(NRP_perc_site_na, pt = 1)

NRP_perc_time <- lm(Perc~SampleEvent, data = NRP_check)
pvalue_table[3,6] <- p_value_ext(NRP_perc_time, pt = 1)

# Res P
mean(soils_PSI$ResP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Acidic")], na.rm = TRUE)*100
sd(soils_PSI$ResP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Acidic")], na.rm = TRUE)/sqrt(length(soils_PSI$ResP[which(soils_PSI$Area=="Acidic")]))*100

mean(soils_PSI$ResP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Non-acidic")], na.rm = TRUE)*100
sd(soils_PSI$ResP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Non-acidic")], na.rm = TRUE)/sqrt(length(soils_PSI$ResP[which(soils_PSI$Area=="Non-acidic")]))*100


resP_check <- soils_PSI %>% dplyr::select(SampleEvent, Area, SiteName, PlotNumber,ResP, PfracTotP)
resP_check$Perc <- resP_check$ResP/resP_check$PfracTotP*100
resP_perc_area <- lm(Perc~Area, data = resP_check)
pvalue_table[4,2] <- p_value_ext(resP_perc_area, pt=1)

resP_perc_site <- lm(Perc~SiteName, data = resP_check)
pvalue_table[4,3] <- p_value_ext(resP_perc_site, pt = 1)

resP_perc_site_a <- lm(Perc~SiteName, data = resP_check[resP_check$Area=="Acidic",])
pvalue_table[4,4] <- p_value_ext(resP_perc_site_a, pt = 1)

resP_perc_site_na <- lm(Perc~SiteName, data = resP_check[resP_check$Area=="Non-acidic",])
pvalue_table[4,5] <- p_value_ext(resP_perc_site_na, pt = 1)

resP_perc_time <- lm(Perc~SampleEvent, data = resP_check)
pvalue_table[4,6] <- p_value_ext(resP_perc_time, pt = 1)


### Reactive Phosphorous fractionation ####
mean(soils_PSI$PfracH2OSRP/soils_PSI$PfracRP)*100
(sd(soils_PSI$PfracH2OSRP/soils_PSI$PfracRP))/sqrt(length(soils_PSI$PfracH2OSRP))*100

## H2OSRP variance ##
H2OSRP_site <- lm(lnH2OSRP~SiteName, data = soils_PSI)
base::summary(H2OSRP_site)
pvalue_table[5,3] <-p_value_ext(H2OSRP_site, pt = 1)

H2OSRP_siteA <- lm(lnH2OSRP~SiteName, data = soils_PSI[soils_PSI$Area=="Acidic",])
base::summary(H2OSRP_siteA)
pvalue_table[5,4] <- p_value_ext(H2OSRP_siteA, pt = 1)

H2OSRP_siteNA <- lm(lnH2OSRP~SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
base::summary(H2OSRP_siteNA)
pvalue_table[5,5] <- p_value_ext(H2OSRP_siteNA, pt = 1)

mean(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Mesic")])
sd(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Mesic")])/sqrt(length(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Mesic")]))
mean(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Dry")])
sd(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Dry")])/sqrt(length(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Dry")]))
mean(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Hydric")])
sd(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Hydric")])/sqrt(length(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Hydric")]))

H2OSRP_area <- lm(lnH2OSRP~Area, data = soils_PSI)
pvalue_table[5,2] <- p_value_ext(H2OSRP_area, pt = 1)

H2OSRP_time <- lm(lnH2OSRP~SampleEvent, data = soils_PSI)
base::summary(H2OSRP_time)
pvalue_table[5,6] <- p_value_ext(H2OSRP_time, pt = 1)

## BDSRP variance #
mean(soils_PSI$PfracBDSRP/soils_PSI$PfracRP)*100
(sd(soils_PSI$PfracBDSRP/soils_PSI$PfracRP))/sqrt(length(soils_PSI$PfracBDSRP))*100

mean(soils_PSI$PfracBDSRP[soils_PSI$Area=="Acidic"]/soils_PSI$PfracRP[soils_PSI$Area=="Acidic"])*100
(sd(soils_PSI$PfracBDSRP[soils_PSI$Area=="Acidic"]/soils_PSI$PfracRP[soils_PSI$Area=="Acidic"]))/sqrt(length(soils_PSI$PfracBDSRP[soils_PSI$Area=="Acidic"]))*100

mean(soils_PSI$PfracBDSRP[soils_PSI$Area=="Non-acidic"]/soils_PSI$PfracRP[soils_PSI$Area=="Non-acidic"])*100
(sd(soils_PSI$PfracBDSRP[soils_PSI$Area=="Non-acidic"]/soils_PSI$PfracRP[soils_PSI$Area=="Non-acidic"]))/sqrt(length(soils_PSI$PfracBDSRP[soils_PSI$Area=="Non-acidic"]))*100

mean(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Mesic")])
sd(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Mesic")])/sqrt(length(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Mesic")]))
mean(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Dry")])
sd(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Dry")])/sqrt(length(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Dry")]))
mean(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Hydric")])
sd(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Hydric")])/sqrt(length(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Hydric")]))


BDrp_site <- lm(PfracBDSRP~SiteName, data = soils_PSI)
base::summary(BDrp_site)
pvalue_table[6,3] <- p_value_ext(BDrp_site, pt = 1)

BDrp_siteA <- lm(PfracBDSRP~SiteName, data = soils_PSI[soils_PSI$Area=="Acidic",])
base::summary(BDrp_siteA)
pvalue_table[6,4] <- p_value_ext(BDrp_siteA, pt = 1)

BDrp_siteNA <- lm(PfracBDSRP~SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
base::summary(BDrp_siteNA)
pvalue_table[6,5] <- p_value_ext(BDrp_siteNA, pt = 1)

BDrp_area <- lm(PfracBDSRP~Area, data = soils_PSI)
base::summary(BDrp_area)
pvalue_table[6,2]<- p_value_ext(BDrp_area, pt = 1)

BDrp_time <- lm(sqrt(PfracBDSRP)~SampleEvent, data = soils_PSI)
base::summary(BDrp_time)
pvalue_table[6,6] <- p_value_ext(BDrp_time, pt = 1)

mean(soils_PSI$PfracBDSRP[soils_PSI$SampleEvent=="Early Season"])
sd(soils_PSI$PfracBDSRP[soils_PSI$SampleEvent=="Early Season"])/sqrt(length(soils_PSI$PfracBDSRP[soils_PSI$SampleEvent=="Early Season"]))
mean(soils_PSI$PfracBDSRP[soils_PSI$SampleEvent=="Mid Season"])
sd(soils_PSI$PfracBDSRP[soils_PSI$SampleEvent=="Mid Season"])/sqrt(length(soils_PSI$PfracBDSRP[soils_PSI$SampleEvent=="Mid Season"]))
mean(soils_PSI$PfracBDSRP[soils_PSI$SampleEvent=="Late Season"])
sd(soils_PSI$PfracBDSRP[soils_PSI$SampleEvent=="Late Season"])/sqrt(length(soils_PSI$PfracBDSRP[soils_PSI$SampleEvent=="Late Season"]))


## NaOHSRP variance by area and site##
AlSRP_area <- lm(PfracNaOHSRP~Area, data = soils_PSI)
anova(AlSRP_area)
shapiro.test(resid(AlSRP_area))
pvalue_table[7,2] <- p_value_ext(AlSRP_area, pt = 1)

AlSRP_treat <- lm(PfracNaOHSRP~treatment, data = soils_PSI)
anova(AlSRP_treat)
shapiro.test(resid(AlSRP_treat))
base::summary(AlSRP_treat)

AlSRP_site <- lm(PfracNaOHSRP~SiteName, data = soils_PSI)
shapiro.test(resid(AlSRP_site))
anova(AlSRP_site)
base::summary(AlSRP_site)
pvalue_table[7,3]<- p_value_ext(AlSRP_site, pt = 1)

AlSRP_siteA <- lm(PfracNaOHSRP~SiteName, data = soils_PSI[soils_PSI$Area=="Acidic",])
shapiro.test(resid(AlSRP_siteA))
anova(AlSRP_siteA)
base::summary(AlSRP_siteA)
pvalue_table[7,4]<- p_value_ext(AlSRP_siteA, pt = 1)

AlSRP_siteNA <- lm(PfracNaOHSRP~SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(AlSRP_siteNA))
anova(AlSRP_siteNA)
base::summary(AlSRP_siteNA)
pvalue_table[7,5]<- p_value_ext(AlSRP_siteNA, pt = 1)

mean(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Acidic")])*100
(sd(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Acidic")]))/sqrt(length(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Acidic")]))*100
mean(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Non-acidic")])*100
(sd(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Non-acidic")]))/sqrt(length(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Non-acidic")]))*100


AlSRP_time <- lm(PfracNaOHSRP~SampleEvent, data = soils_PSI)
anova(AlSRP_time)
shapiro.test(resid(AlSRP_time))
base::summary(AlSRP_time)
pvalue_table[7,6]<- p_value_ext(AlSRP_time, pt = 1)

## HClP variance by area and site##
CaP_area <- lm(lnHClP~Area, data = soils_PSI)
anova(CaP_area)
shapiro.test(resid(CaP_area))
pvalue_table[8,2] <- p_value_ext(CaP_area, pt = 1)

CaP_treat <- lm(lnHClP~treatment, data = soils_PSI)
anova(CaP_treat)
shapiro.test(resid(CaP_treat))
base::summary(CaP_treat)

CaP_site <- lm(lnHClP~SiteName, data = soils_PSI)
anova(CaP_site)
shapiro.test(resid(CaP_site))
base::summary(CaP_site)
pvalue_table[8,3] <- p_value_ext(CaP_site, pt = 1)

CaP_siteA <- lm(lnHClP~SiteName, data = soils_PSI[soils_PSI$Area=="Acidic",])
anova(CaP_siteA)
shapiro.test(resid(CaP_siteA))
base::summary(CaP_siteA)
pvalue_table[8,4] <- p_value_ext(CaP_siteA, pt = 1)

CaP_siteNA <- lm(lnHClP~SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
anova(CaP_siteNA)
shapiro.test(resid(CaP_siteNA))
base::summary(CaP_siteNA)
pvalue_table[8,5] <- p_value_ext(CaP_siteNA, pt = 1)

mean(soils_PSI$PfracHClP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Acidic")])*100
(sd(soils_PSI$PfracHClP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Acidic")]))/sqrt(length(soils_PSI$PfracHClP[which(soils_PSI$Area=="Acidic")]))*100
mean(soils_PSI$PfracHClP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Non-acidic")])*100
(sd(soils_PSI$PfracHClP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Non-acidic")]))/sqrt(length(soils_PSI$PfracHClP[which(soils_PSI$Area=="Non-acidic")]))*100

CaP_time <- lm(lnHClP~SampleEvent, data = soils_PSI)
anova(CaP_time)
shapiro.test(resid(CaP_time))
base::summary(CaP_time)
pvalue_table[8,6] <- p_value_ext(CaP_time, pt = 1)

### Iron fractionation ####
## Total Fe ###
Fe_area <- lm(sqrt(PfracTotFe)~Area, data= soils_PSI)
base::summary(Fe_area)
anova(Fe_area)
shapiro.test(resid(Fe_area))
pvalue_table[9,2] <- p_value_ext(Fe_area, pt = 1)

Fe_site <- lm(sqrt(PfracTotFe)~SiteName, data= soils_PSI)
base::summary(Fe_site)
anova(Fe_site)
shapiro.test(resid(Fe_site))
pvalue_table[9,3] <- p_value_ext(Fe_site, pt = 1)

Fe_time <- lm(sqrt(PfracTotFe)~SampleEvent, data= soils_PSI)
base::summary(Fe_time)
anova(Fe_time)
shapiro.test(resid(Fe_time))
pvalue_table[9,6] <- p_value_ext(Fe_time, pt = 1)

Fe_a_site <- lm(sqrt(PfracTotFe)~SiteName, data= soils_PSI[which(soils_PSI$Area=="Acidic"),])
base::summary(Fe_a_site)
anova(Fe_a_site)
shapiro.test(resid(Fe_a_site))
pvalue_table[9,4] <- p_value_ext(Fe_a_site, pt = 1)

Fe_na_site <- lm(sqrt(PfracTotFe)~SiteName, data= soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
base::summary(Fe_na_site)
anova(Fe_na_site)
shapiro.test(resid(Fe_na_site))
pvalue_table[9,5] <- p_value_ext(Fe_na_site, pt = 1)

mean(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Dry")])
sd(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Dry")])/sqrt(length(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Dry")]))

mean(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Mesic")])
sd(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Mesic")])/sqrt(length(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Mesic")]))

mean(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Hydric")])
sd(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Hydric")])/sqrt(length(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Hydric")]))

mean(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Dry")])
sd(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Dry")])/sqrt(length(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Dry")]))

mean(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Mesic")])
sd(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Mesic")])/sqrt(length(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Mesic")]))

mean(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Hydric")])
sd(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Hydric")])/sqrt(length(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Hydric")]))

# organic-bound Fe
shapiro.test(sqrt(soils_PSI$OrgFe))
OFe_area <- lm(sqrt(OrgFe)~Area, data = soils_PSI)
shapiro.test(resid(OFe_area))
base::summary(OFe_area)
pvalue_table[10,2] <- p_value_ext(OFe_area, pt = 1)

OFe_a_site <- lm(sqrt(OrgFe)~SiteName, data = soils_PSI[which(soils_PSI$Area=="Acidic"),])
shapiro.test(resid(OFe_a_site))
base::summary(OFe_a_site)
pvalue_table[10,4] <- p_value_ext(OFe_a_site, pt = 1)

OFe_na_site <- lm(sqrt(OrgFe)~SiteName, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
shapiro.test(resid(OFe_na_site))
base::summary(OFe_na_site)
pvalue_table[10,5] <- p_value_ext(OFe_na_site, pt = 1)

OFe_site <- lm(sqrt(OrgFe)~SiteName, data = soils_PSI)
shapiro.test(resid(OFe_site))
base::summary(OFe_site)
pvalue_table[10,3] <- p_value_ext(OFe_site, pt = 1)

OFe_time <- lm(sqrt(OrgFe)~SampleEvent, data = soils_PSI)
shapiro.test(resid(OFe_time))
base::summary(OFe_time)
pvalue_table[10,6] <- p_value_ext(OFe_time, pt = 1)

mean(soils_PSI$OrgFe[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Acidic")])*100
sd(soils_PSI$OrgFe[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Acidic")])/sqrt(length(soils_PSI$OrgFe[which(soils_PSI$Area=="Acidic")]))*100

mean(soils_PSI$OrgFe[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Non-acidic")])*100
(sd(soils_PSI$OrgFe[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Non-acidic")])/sqrt(length(soils_PSI$OrgFe[which(soils_PSI$Area=="Non-acidic")])))*100

mean(soils_PSI$OrgFe[soils_PSI$Area=="Acidic"])
sd(soils_PSI$OrgFe[soils_PSI$Area=="Acidic"])/sqrt(length(soils_PSI$OrgFe[soils_PSI$Area=="Acidic"]))

mean(soils_PSI$OrgFe[soils_PSI$Area=="Non-acidic"])
sd(soils_PSI$OrgFe[soils_PSI$Area=="Non-acidic"])/sqrt(length(soils_PSI$OrgFe[soils_PSI$Area=="Non-acidic"]))

# non-crystalline iron
NC_Fe_area <- lm(sq_Fe~Area, data = soils_PSI)
shapiro.test(resid(NC_Fe_area))
base::summary(NC_Fe_area)
pvalue_table[11,2] <- p_value_ext(NC_Fe_area, pt = 1)

mean(soils_PSI$NonCFe[soils_PSI$Area=="Non-acidic"])
sd(soils_PSI$NonCFe[soils_PSI$Area=="Non-acidic"])/sqrt(length(soils_PSI$NonCFe[soils_PSI$Area=="Non-acidic"]))
mean(soils_PSI$NonCFe[soils_PSI$Area=="Acidic"])
sd(soils_PSI$NonCFe[soils_PSI$Area=="Acidic"])/sqrt(length(soils_PSI$NonCFe[soils_PSI$Area=="Acidic"]))

mean(soils_PSI$NonCFe[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Acidic")])*100
(sd(soils_PSI$NonCFe[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Acidic")])/sqrt(length(soils_PSI$NonCFe[which(soils_PSI$Area=="Acidic")])))*100

mean(soils_PSI$FefracDHFe[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Acidic")])*100
(sd(soils_PSI$FefracDHFe[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Acidic")])/sqrt(length(soils_PSI$FefracDHFe[which(soils_PSI$Area=="Acidic")])))*100

mean(soils_PSI$NonCFe[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Non-acidic")])*100
(sd(soils_PSI$NonCFe[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Non-acidic")])/sqrt(length(soils_PSI$NonCFe[which(soils_PSI$Area=="Non-acidic")])))*100

mean(soils_PSI$FefracDHFe[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Non-acidic")])*100
(sd(soils_PSI$FefracDHFe[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Non-acidic")])/sqrt(length(soils_PSI$FefracDHFe[which(soils_PSI$Area=="Non-acidic")])))*100

mean(soils_PSI$NonCFe[which(soils_PSI$treatment=="Acidic Hydric")])
sd(soils_PSI$NonCFe[which(soils_PSI$treatment=="Acidic Hydric")])/sqrt(length(soils_PSI$NonCFe[which(soils_PSI$treatment=="Acidic Hydric")]))
mean(soils_PSI$NonCFe[which(soils_PSI$treatment=="Acidic Dry")])
sd(soils_PSI$NonCFe[which(soils_PSI$treatment=="Acidic Dry")])/sqrt(length(soils_PSI$NonCFe[which(soils_PSI$treatment=="Acidic Dry")]))
mean(soils_PSI$NonCFe[which(soils_PSI$treatment=="Acidic Mesic")])
sd(soils_PSI$NonCFe[which(soils_PSI$treatment=="Acidic Mesic")])/sqrt(length(soils_PSI$NonCFe[which(soils_PSI$treatment=="Acidic Mesic")]))

NC_Fe_a_site <- lm(sq_Fe~SiteName, data = soils_PSI[which(soils_PSI$Area=="Acidic"),])
shapiro.test(resid(NC_Fe_a_site))
base::summary(NC_Fe_a_site)
pvalue_table[11,4] <- p_value_ext(NC_Fe_a_site, pt = 1)

NC_Fe_na_site <- lm(sq_Fe~SiteName, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
shapiro.test(resid(NC_Fe_na_site))
base::summary(NC_Fe_na_site)
pvalue_table[11,5] <- p_value_ext(NC_Fe_na_site, pt = 1)

NC_Fe_site <- lm(sq_Fe~SiteName, data = soils_PSI)
shapiro.test(resid(NC_Fe_site))
base::summary(NC_Fe_site)
pvalue_table[11,3] <- p_value_ext(NC_Fe_site, pt = 1)

NC_Fe_time <- lm(sq_Fe~SampleEvent, data = soils_PSI)
shapiro.test(resid(NC_Fe_time))
base::summary(NC_Fe_time)
pvalue_table[11,6] <- p_value_ext(NC_Fe_time, pt = 1)

## organic bound Fe and non-crystalline iron ##
mean(soils_PSI$OrgFe/soils_PSI$NonCFe)*100
sd(soils_PSI$OrgFe/soils_PSI$NonCFe)/sqrt(length(soils_PSI$NonCFe))*100

mean(soils_PSI$OrgFe[soils_PSI$SiteName=="Mesic"]/soils_PSI$NonCFe[soils_PSI$SiteName=="Mesic"])*100
sd(soils_PSI$OrgFe[soils_PSI$SiteName=="Mesic"]/soils_PSI$NonCFe[soils_PSI$SiteName=="Mesic"])/sqrt(length(soils_PSI$NonCFe))*100

## iron oxides ##
BDFe_area <- lm(log(PfracBDFe)~Area, data = soils_PSI)
shapiro.test(resid(BDFe_area))
anova(BDFe_area)
base::summary(BDFe_area)
pvalue_table[13,2] <- p_value_ext(BDFe_area, pt = 1)

BDFe_siteA <- lm(log(PfracBDFe)~SiteName, data = soils_PSI[which(soils_PSI$Area=="Acidic"),])
shapiro.test(resid(BDFe_siteA))
base::summary(BDFe_siteA)
pvalue_table[13,4] <- p_value_ext(BDFe_siteA, pt = 1)

BDFe_siteNA <- lm(log(PfracBDFe)~SiteName, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
shapiro.test(resid(BDFe_siteNA))
base::summary(BDFe_siteNA)
pvalue_table[13,5] <- p_value_ext(BDFe_siteNA, pt = 1)

BDFe_site <- lm(log(PfracBDFe)~SiteName, data = soils_PSI)
shapiro.test(resid(BDFe_site))
base::summary(BDFe_site)
pvalue_table[13,3] <- p_value_ext(BDFe_site, pt = 1)

BDFe_time <- lm(log(PfracBDFe)~SampleEvent, data = soils_PSI)
shapiro.test(resid(BDFe_time))
base::summary(BDFe_time)
pvalue_table[13,6] <- p_value_ext(BDFe_time, pt = 1)


# iron oxides and total extracted iron
BD_NC_cor <- cor.test(log(soils_PSI$PfracBDFe), log(soils_PSI$NonCFe+soils_PSI$FefracDHFe), method = "pearson")
BD_NC <- lm(log(PfracBDFe)~ log(NonCFe+FefracDHFe), data = soils_PSI)
shapiro.test(resid(BD_NC))
pvalue_table_pt2[2,2] <-p_value_ext(BD_NC,pt=2,cor = BD_NC_cor)

BD_NC_site <- lm(log(PfracBDFe)~ log(NonCFe+FefracDHFe)*SiteName, data = soils_PSI)
shapiro.test(resid(BD_NC_site))
pvalue_table_pt2[2,3] <-p_value_ext(BD_NC_site,pt=1)

## crystalline iron ##
DHFe_Area <- lm(log(FefracDHFe)~Area, data = soils_PSI)
shapiro.test(resid(DHFe_Area))
base::summary(DHFe_Area)
pvalue_table[12,2] <- p_value_ext(DHFe_Area, pt = 1)

DHFe_siteA <- lm(log(FefracDHFe)~SiteName, data = soils_PSI[soils_PSI$Area=="Acidic",])
shapiro.test(resid(DHFe_siteA))
base::summary(DHFe_siteA)
pvalue_table[12,4] <- p_value_ext(DHFe_siteA, pt = 1)

DHFe_siteNA <- lm(log(FefracDHFe)~SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(DHFe_siteNA))
base::summary(DHFe_siteNA)
pvalue_table[12,5] <- p_value_ext(DHFe_siteNA, pt = 1)

DHFe_site <- lm(log(FefracDHFe)~SiteName, data = soils_PSI)
shapiro.test(resid(DHFe_site))
base::summary(DHFe_site)
pvalue_table[12,3] <- p_value_ext(DHFe_site, pt = 1)

DHFe_time <- lm(log(FefracDHFe)~SampleEvent, data = soils_PSI)
shapiro.test(resid(DHFe_time))
base::summary(DHFe_time)
pvalue_table[12,6] <- p_value_ext(DHFe_time, pt = 1)

### reactive phosphorus pools ####
## BDSRP covariance with H2OSRP ##
SRP_H2O <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI)
anova(SRP_H2O)
base::summary(SRP_H2O)
shapiro.test(resid(SRP_H2O))
SRP_H2O_cor <- cor.test(soils_PSI$lnH2OSRP,soils_PSI$PfracBDSRP, method = "pearson")
pvalue_table_pt3[1,2] <- p_value_ext(SRP_H2O, pt = 2,cor = SRP_H2O_cor)
pvalue_table_pt2[1,2] <- p_value_ext(SRP_H2O, pt = 2,cor = SRP_H2O_cor)

ancova(lnH2OSRP~PfracBDSRP*SiteName, data = soils_PSI)
SRP_H2O_site <- lm(lnH2OSRP~PfracBDSRP*SiteName, data = soils_PSI)
base::summary(SRP_H2O_site)
pvalue_table_pt2[1,3] <- p_value_ext(SRP_H2O_site, pt = 1)

SRP_H2O_D <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI[which(soils_PSI$SiteName=="Dry"),])
base::summary(SRP_H2O_D)
coef(SRP_H2O_D)
SRP_H2O_D_cor <- cor.test(soils_PSI$lnH2OSRP[which(soils_PSI$SiteName=="Dry")], soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Dry")])
pvalue_table_pt3[1,3] <- p_value_ext(SRP_H2O_D, pt = 2,cor = SRP_H2O_D_cor)

SRP_H2O_M <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI[which(soils_PSI$SiteName=="Mesic"),])
base::summary(SRP_H2O_M)
coef(SRP_H2O_M)
SRP_H2O_M_cor <- cor.test(soils_PSI$lnH2OSRP[which(soils_PSI$SiteName=="Mesic")], soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Mesic")])
pvalue_table_pt3[1,4] <- p_value_ext(SRP_H2O_M, pt = 2, cor = SRP_H2O_M_cor)

SRP_H2O_H <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI[which(soils_PSI$SiteName=="Hydric"),])
base::summary(SRP_H2O_H)
coef(SRP_H2O_H)
SRP_H2O_H_cor <- cor.test(soils_PSI$lnH2OSRP[which(soils_PSI$SiteName=="Hydric")], soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Hydric")])
pvalue_table_pt3[1,5] <- p_value_ext(SRP_H2O_H, pt = 2, cor = SRP_H2O_H_cor)

## SRP covariance between H2O, NaOH, and HCl ##
NaOH_H2O <- lm(lnH2OSRP~PfracNaOHSRP, data = soils_PSI)
NaOH_H2O_cor <- cor.test(soils_PSI$lnH2OSRP,soils_PSI$PfracNaOHSRP, method = "pearson")
pvalue_table_pt3[2,2] <- p_value_ext(NaOH_H2O, pt = 2, NaOH_H2O_cor)

ancova(lnH2OSRP~PfracNaOHSRP*SiteName, data = soils_PSI)

NaOH_H2O_D <- lm(lnH2OSRP~PfracNaOHSRP, data = soils_PSI[which(soils_PSI$SiteName=="Dry"),])
base::summary(NaOH_H2O_D)
NaOH_H2O_D_cor <- cor.test(soils_PSI$lnH2OSRP[which(soils_PSI$SiteName=="Dry")], soils_PSI$PfracNaOHSRP[which(soils_PSI$SiteName=="Dry")], )
pvalue_table_pt3[2,3] <- p_value_ext(NaOH_H2O_D, pt = 2, cor = NaOH_H2O_D_cor)

NaOH_H2O_M <- lm(lnH2OSRP~PfracNaOHSRP, data = soils_PSI[which(soils_PSI$SiteName=="Mesic"),])
base::summary(NaOH_H2O_M)
NaOH_H2O_M_cor <- cor.test(soils_PSI$PfracNaOHSRP[which(soils_PSI$SiteName == "Mesic")], soils_PSI$lnH2OSRP[which(soils_PSI$SiteName == "Mesic")], method = "pearson")
pvalue_table_pt3[2,4] <- p_value_ext(NaOH_H2O_M, pt = 2, cor = NaOH_H2O_M_cor)

NaOH_H2O_H <- lm(lnH2OSRP~PfracNaOHSRP, data = soils_PSI[which(soils_PSI$SiteName=="Hydric"),])
base::summary(NaOH_H2O_H)
NaOH_H2O_H_cor <- cor.test(soils_PSI$PfracNaOHSRP[which(soils_PSI$SiteName == "Hydric")], soils_PSI$lnH2OSRP[which(soils_PSI$SiteName == "Hydric")], method = "pearson")
pvalue_table_pt3[2,5] <- p_value_ext(NaOH_H2O_H, pt = 2, cor = NaOH_H2O_H_cor)


HCl_H2O <-lm(lnH2OSRP~PfracHClP, data = soils_PSI)
HCl_H2O_cor <- cor.test(soils_PSI$lnH2OSRP, soils_PSI$PfracHClP, method = "pearson")
pvalue_table_pt3[3,2] <- p_value_ext(HCl_H2O, pt = 2, cor = HCl_H2O_cor)

HCl_H2O_D <-lm(lnH2OSRP~PfracHClP, data = soils_PSI[soils_PSI$SiteName=="Dry",])
HCl_H2O_cor_D <- cor.test(soils_PSI$lnH2OSRP[soils_PSI$SiteName=="Dry"], soils_PSI$PfracHClP[soils_PSI$SiteName=="Dry"], method = "pearson")
pvalue_table_pt3[3,3] <- p_value_ext(HCl_H2O_D, pt = 2, cor = HCl_H2O_cor_D)

HCl_H2O_M <-lm(lnH2OSRP~PfracHClP, data = soils_PSI[soils_PSI$SiteName=="Mesic",])
HCl_H2O_cor_M <- cor.test(soils_PSI$lnH2OSRP[soils_PSI$SiteName=="Mesic"], soils_PSI$PfracHClP[soils_PSI$SiteName=="Mesic"], method = "pearson")
pvalue_table_pt3[3,4] <- p_value_ext(HCl_H2O_M, pt = 2, cor = HCl_H2O_cor_M)

HCl_H2O_H <-lm(lnH2OSRP~PfracHClP, data = soils_PSI[soils_PSI$SiteName=="Hydric",])
HCl_H2O_cor_H <- cor.test(soils_PSI$lnH2OSRP[soils_PSI$SiteName=="Hydric"], soils_PSI$PfracHClP[soils_PSI$SiteName=="Hydric"], method = "pearson")
pvalue_table_pt3[3,5] <- p_value_ext(HCl_H2O_H, pt = 2, cor = HCl_H2O_cor_H)

### H2OSRP and minerals ####
H2Orp_OFe <- lm(PfracH2OSRP~OrgFe, data = soils_PSI)
cor.test(soils_PSI$PfracH2OSRP, soils_PSI$OrgFe, method = "spearman")
pvalue_table_pt2[3,2] <- p_value_ext(H2Orp_OFe, pt = 3,H2Orp_OFe_cor)

H2Orp_OFe_site <- lm(PfracH2OSRP~OrgFe*SiteName, data = soils_PSI)
shapiro.test(resid(H2Orp_OFe_site))
anova(H2Orp_OFe_site)
pvalue_table_pt2[3,3] <- p_value_ext(H2Orp_OFe_site, pt = 1)

## H2OSRP and iron oxides##
H2Orp_BDFe <- lm(PfracH2OSRP~PfracBDFe, data = soils_PSI)
shapiro.test(resid(H2Orp_BDFe))
anova(H2Orp_BDFe)
H2Orp_BDFe_cor <- cor.test(soils_PSI$PfracH2OSRP, soils_PSI$PfracBDFe, method = "spearman")
# r = -0.27, p = 0.05
pvalue_table_pt2[4,2] <- p_value_ext(H2Orp_BDFe, pt = 3, cor = H2Orp_BDFe_cor)

H2Orp_BDFe_site <- lm(PfracH2OSRP~PfracBDFe*SiteName, data = soils_PSI)
shapiro.test(resid(H2Orp_BDFe_site))
anova(H2Orp_BDFe_site)
pvalue_table_pt2[4,3] <- p_value_ext(H2Orp_BDFe_site, pt = 1)

# H2OSRP and non-crystalline
H2Orp_NCFe <- lm(PfracH2OSRP~NonCFe, data = soils_PSI)
shapiro.test(resid(H2Orp_NCFe))
anova(H2Orp_NCFe)
H2Orp_NCFe_cor <- cor.test(soils_PSI$PfracH2OSRP, soils_PSI$NonCFe, method = "spearman")
# r = -0.21, p = 0.15
pvalue_table_pt2[5,2] <- p_value_ext(H2Orp_NCFe, pt = 3, cor = H2Orp_NCFe_cor)

H2Orp_NCFe_site <- lm(PfracH2OSRP~NonCFe*SiteName, data = soils_PSI)
shapiro.test(resid(H2Orp_NCFe_site))
anova(H2Orp_NCFe_site)
pvalue_table_pt2[5,3] <- p_value_ext(H2Orp_NCFe_site, pt = 1)

# H2OSRP and crystalline iron
H2Orp_DHFe <- lm(PfracH2OSRP~FefracDHFe, data = soils_PSI)
shapiro.test(resid(H2Orp_DHFe))
anova(H2Orp_DHFe)
H2Orp_DHFe_cor <- cor.test(soils_PSI$PfracH2OSRP, soils_PSI$FefracDHFe, method = "spearman")
# r = -0.1, p = 0.46
pvalue_table_pt2[6,2] <- p_value_ext(H2Orp_DHFe, pt = 3, cor = H2Orp_DHFe_cor)

H2Orp_DHFe_site <- lm(PfracH2OSRP~FefracDHFe*SiteName, data = soils_PSI)
shapiro.test(resid(H2Orp_DHFe_site))
anova(H2Orp_DHFe_site)
pvalue_table_pt2[6,3] <- p_value_ext(H2Orp_DHFe_site, pt = 1)

## H2O~SRP vs Ca and Al ##
H2O_Al_a <- lm(PfracH2OSRP~PfracNaOHAl, data = soils_PSI[soils_PSI$Area=="Acidic",])
shapiro.test(resid(H2O_Al_a))
anova(H2O_Al_a)
base::summary(H2O_Al_a)
H2O_Al_a_cor <- cor.test(soils_PSI$PfracH2OSRP[soils_PSI$Area=="Acidic"], soils_PSI$PfracNaOHAl[soils_PSI$Area=="Acidic"], method = "pearson")
pvalue_table_pt2[7,2] <- p_value_ext(H2O_Al_a, pt = 2, H2O_Al_a_cor)

H2O_Al_a_site <- lm(PfracH2OSRP~PfracNaOHAl*SiteName, data = soils_PSI[soils_PSI$Area=="Acidic",])
shapiro.test(resid(H2O_Al_a_site))
anova(H2O_Al_a_site)
pvalue_table_pt2[7,3] <- p_value_ext(H2O_Al_a_site, pt = 1)

H2O_Al_na <- lm(PfracH2OSRP~PfracNaOHAl, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(H2O_Al_na))
anova(H2O_Al_na)
base::summary(H2O_Al_na)
H2O_Al_na_cor <- cor.test(soils_PSI$PfracH2OSRP[soils_PSI$Area=="Non-acidic"], soils_PSI$PfracNaOHAl[soils_PSI$Area=="Non-acidic"], method = "pearson")
pvalue_table_pt2[8,2] <- p_value_ext(H2O_Al_na, pt = 2, H2O_Al_na_cor)

H2O_Al_na_site <- lm(PfracH2OSRP~PfracNaOHAl*SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(H2O_Al_na_site))
anova(H2O_Al_na_site)
pvalue_table_pt2[8,3] <- p_value_ext(H2O_Al_na_site, pt = 1)


H2O_Ca_na <- lm(lnH2OSRP~PfracHClCa, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(H2O_Ca_na))
ancova(PfracH2OSRP~PfracHClCa*Area, data = soils_PSI)
base::summary(H2O_Ca_na)
H2O_Ca_cor_na <- cor.test(soils_PSI$lnH2OSRP[soils_PSI$Area=="Non-acidic"], soils_PSI$PfracHClCa[soils_PSI$Area=="Non-acidic"])
pvalue_table_pt3[7,2] <- p_value_ext(H2O_Ca_na, pt = 2, cor = H2O_Ca_cor_na)

ancova(PfracH2OSRP~PfracHClCa*SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])

H2OSRP_Ca_NAD_cor <- cor.test(soils_PSI$PfracH2OSRP[which(soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Dry")],soils_PSI$PfracHClCa[which(soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Dry")], method = "pearson")
H2OSRP_Ca_NAD <- lm(PfracH2OSRP~PfracHClCa, data = soils_PSI[soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Dry",])
pvalue_table_pt3[7,3] <- p_value_ext(H2OSRP_Ca_NAD, pt = 2, H2OSRP_Ca_NAD_cor)

H2OSRP_Ca_NAM_cor <- cor.test(soils_PSI$PfracH2OSRP[which(soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Mesic")],soils_PSI$PfracHClCa[which(soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Mesic")], method = "pearson")
H2OSRP_Ca_NAM <- lm(PfracH2OSRP~PfracHClCa, data = soils_PSI[soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Mesic",])
pvalue_table_pt3[7,4] <- p_value_ext(H2OSRP_Ca_NAM, pt = 2, H2OSRP_Ca_NAM_cor)

H2OSRP_Ca_NAH_cor <- cor.test(soils_PSI$PfracH2OSRP[which(soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Hydric")],soils_PSI$PfracHClCa[which(soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Hydric")], method = "pearson")
H2OSRP_Ca_NAH <- lm(PfracH2OSRP~PfracHClCa, data = soils_PSI[soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Hydric",])
pvalue_table_pt3[7,5] <- p_value_ext(H2OSRP_Ca_NAH, pt = 2, H2OSRP_Ca_NAH_cor)

### BD~SRP and iron ####
## BD~SRP and Organic bound Fe ##
OFe_BDSRP <- lm(PfracBDSRP~OrgFe, data = soils_PSI)
shapiro.test(resid(OFe_BDSRP))
anova(OFe_BDSRP)         

OFe_BDSRP_site <- lm(PfracBDSRP~OrgFe*Area, data = soils_PSI)
shapiro.test(resid(OFe_BDSRP_site))
anova(OFe_BDSRP_site)
base::summary(OFe_BDSRP_site)

OFe_BDSRP_a_cor <- cor.test(soils_PSI$PfracBDSRP[soils_PSI$Area =="Acidic"],soils_PSI$OrgFe[soils_PSI$Area=="Acidic"])
OFe_BDSRP_a <- lm(PfracBDSRP~OrgFe, data = soils_PSI)
pvalue_table_pt3[5,2] <- p_value_ext(OFe_BDSRP_a, pt=2,OFe_BDSRP_a_cor)

OFe_BDSRP_AD_cor <- cor.test(soils_PSI$PfracBDSRP[soils_PSI$treatment =="Acidic Dry"],soils_PSI$OrgFe[soils_PSI$treatment=="Acidic Dry"])
pvalue_table_pt3[5,3] <- p_value_ext(pt =2,cor = OFe_BDSRP_AD_cor)

OFe_BDSRP_AM_cor <- cor.test(soils_PSI$PfracBDSRP[soils_PSI$treatment =="Acidic Mesic"],soils_PSI$OrgFe[soils_PSI$treatment=="Acidic Mesic"])
pvalue_table_pt3[5,4] <- p_value_ext(pt =2,cor = OFe_BDSRP_AM_cor)

OFe_BDSRP_AH_cor <- cor.test(soils_PSI$PfracBDSRP[soils_PSI$treatment =="Acidic Hydric"],soils_PSI$OrgFe[soils_PSI$treatment=="Acidic Hydric"])
pvalue_table_pt3[5,5] <- p_value_ext(pt =2,cor = OFe_BDSRP_AH_cor)


## BD~SRP and non-crystalline iron ##
NCFe_BDSRP <- lm(PfracBDSRP~NonCFe, data = soils_PSI)
shapiro.test(resid(NCFe_BDSRP))
anova(NCFe_BDSRP)
NCFe_BDSRP_cor <- cor.test(soils_PSI$PfracBDSRP,soils_PSI$NonCFe)
pvalue_table_pt3[4,2] <- p_value_ext(NCFe_BDSRP, pt = 2, NCFe_BDSRP_cor)

NCFe_BDSRP_site <- lm(PfracBDSRP~NonCFe*SiteName, data = soils_PSI)
shapiro.test(resid(NCFe_BDSRP_site))
anova(NCFe_BDSRP_site)
base::summary(NCFe_BDSRP_site)

BDSRP_NCFe_D <- lm(PfracBDSRP~log(NonCFe), data = soils_PSI[soils_PSI$SiteName=="Dry",])
anova(BDSRP_NCFe_D)
BDSRP_NCFe_D_cor <- cor.test(soils_PSI$PfracBDSRP[soils_PSI$SiteName=="Dry"],log(soils_PSI$NonCFe[soils_PSI$SiteName=="Dry"]))
pvalue_table_pt3[4,3] <- p_value_ext(BDSRP_NCFe_D, pt = 2, cor = BDSRP_NCFe_D_cor)

BDSRP_NCFe_M <- lm(PfracBDSRP~log(NonCFe), data = soils_PSI[soils_PSI$SiteName=="Mesic",])
anova(BDSRP_NCFe_M)
BDSRP_NCFe_M_cor <- cor.test(soils_PSI$PfracBDSRP[soils_PSI$SiteName=="Mesic"],log(soils_PSI$NonCFe[soils_PSI$SiteName=="Mesic"]))
pvalue_table_pt3[4,4] <- p_value_ext(BDSRP_NCFe_M, pt = 2, cor = BDSRP_NCFe_M_cor)

BDSRP_NCFe_H <- lm(PfracBDSRP~NonCFe, data = soils_PSI[soils_PSI$SiteName=="Hydric",])
anova(BDSRP_NCFe_H)
BDSRP_NCFe_H_cor <-cor.test(soils_PSI$PfracBDSRP[soils_PSI$SiteName=="Hydric"],log(soils_PSI$NonCFe[soils_PSI$SiteName=="Hydric"]))
pvalue_table_pt3[4,5] <- p_value_ext(BDSRP_NCFe_H, pt = 2, cor = BDSRP_NCFe_H_cor)

## BDSRP~crystalline iron ##
BDSRP_DH_a <- lm(PfracBDSRP~FefracDHFe, data = soils_PSI[soils_PSI$Area=="Acidic",])
BDSRP_DH_a_cor <- cor.test(soils_PSI$PfracBDSRP[soils_PSI$Area=="Acidic"], soils_PSI$FefracDHFe[soils_PSI$Area=="Acidic"])
pvalue_table_pt3[6,2] <- p_value_ext(BDSRP_DH_a, pt = 2, cor = BDSRP_DH_a_cor)

BDSRP_DHFe_AD_cor <- cor.test(soils_PSI[which(soils_PSI$treatment=="Acidic Dry"),]$PfracBDSRP, soils_PSI[which(soils_PSI$treatment=="Acidic Dry"),]$FefracDHFe, method = "pearson")
BDSRP_DHFe_AD <- lm(PfracBDSRP~FefracDHFe, data = soils_PSI[soils_PSI$treatment=="Acidic Dry",])
pvalue_table_pt3[6,3] <- p_value_ext(BDSRP_DHFe_AD, pt = 2, cor = BDSRP_DHFe_AD_cor)

BDSRP_DHFe_AM_cor <- cor.test(soils_PSI[which(soils_PSI$treatment=="Acidic Mesic"),]$PfracBDSRP, soils_PSI[which(soils_PSI$treatment=="Acidic Mesic"),]$FefracDHFe, method = "pearson")
BDSRP_DHFe_AM <- lm(PfracBDSRP~FefracDHFe, data = soils_PSI[soils_PSI$treatment=="Acidic Mesic",])
pvalue_table_pt3[6,4] <- p_value_ext(BDSRP_DHFe_AM, pt = 2, cor = BDSRP_DHFe_AM_cor)

BDSRP_DHFe_AH_cor <- cor.test(soils_PSI[which(soils_PSI$treatment=="Acidic Hydric"),]$PfracBDSRP, soils_PSI[which(soils_PSI$treatment=="Acidic Hydric"),]$FefracDHFe, method = "pearson")
BDSRP_DHFe_AH <- lm(PfracBDSRP~FefracDHFe, data = soils_PSI[soils_PSI$treatment=="Acidic Hydric",])
pvalue_table_pt3[6,5] <- p_value_ext(BDSRP_DHFe_AH, pt = 2, cor = BDSRP_DHFe_AH_cor)

## Aluminum and Calcium  ####
Al_area <- lm(PfracNaOHAl~Area, soils_PSI)
shapiro.test(resid(Al_area))
anova(Al_area)
pvalue_table[14,2] <- p_value_ext(Al_area, pt =1)

Al_siteA <- lm(PfracNaOHAl~SiteName, soils_PSI[soils_PSI$Area=="Acidic",])
shapiro.test(resid(Al_siteA))
anova(Al_siteA)
pvalue_table[14,4] <- p_value_ext(Al_siteA, pt = 1)

Al_siteNA <- lm(PfracNaOHAl~SiteName, soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(Al_siteNA))
anova(Al_siteNA)
pvalue_table[14,5] <- p_value_ext(Al_siteNA, pt =1)

Al_site <- lm(PfracNaOHAl~SiteName, soils_PSI)
shapiro.test(resid(Al_site))
anova(Al_site)
pvalue_table[14,3] <- p_value_ext(Al_site, pt =1)

Al_time <- lm(PfracNaOHAl~SampleEvent, soils_PSI)
shapiro.test(resid(Al_time))
anova(Al_time)
pvalue_table[14,6] <- p_value_ext(Al_time, pt =1)

range(soils_PSI$PfracNaOHAl[which(soils_PSI$Area=="Acidic")])
range(soils_PSI$PfracNaOHAl[which(soils_PSI$Area=="Non-acidic")])

TotAl_area <- lm(log(PfracTotAl)~Area, data = soils_PSI)
shapiro.test(resid(TotAl_area))
anova(TotAl_area)
pvalue_table[16,2] <- p_value_ext(TotAl_area, pt = 1)

TotAl_site <- lm(log(PfracTotAl)~SiteName, data = soils_PSI)
shapiro.test(resid(TotAl_site))
anova(TotAl_site)
pvalue_table[16,3] <- p_value_ext(TotAl_site, pt = 1)

TotAl_site_a <- lm(log(PfracTotAl)~SiteName, data = soils_PSI[soils_PSI$Area=="Acidic",])
shapiro.test(resid(TotAl_site_a))
anova(TotAl_site_a)
pvalue_table[16,4] <- p_value_ext(TotAl_site_a, pt = 1)

TotAl_site_na <- lm(log(PfracTotAl)~SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(TotAl_site_na))
anova(TotAl_site_na)
pvalue_table[16,5] <- p_value_ext(TotAl_site_na, pt = 1)

TotAl_time <- lm(log(PfracTotAl)~SampleEvent, data = soils_PSI)
shapiro.test(resid(TotAl_time))
anova(TotAl_time)
pvalue_table[16,6] <- p_value_ext(TotAl_time, pt = 1)

mean(soils_PSI$PfracTotAl[soils_PSI$Area=="Acidic"])
sd(soils_PSI$PfracTotAl[soils_PSI$Area=="Acidic"])/sqrt(length(soils_PSI$PfracTotAl[soils_PSI$Area=="Acidic"]))

mean(soils_PSI$PfracTotAl[soils_PSI$Area=="Non-acidic"])
sd(soils_PSI$PfracTotAl[soils_PSI$Area=="Non-acidic"])/sqrt(length(soils_PSI$PfracTotAl[soils_PSI$Area=="Non-acidic"]))

mean(soils_PSI$PfracNaOHAl/soils_PSI$PfracTotAl)*100
sd(soils_PSI$PfracNaOHAl/soils_PSI$PfracTotAl)/sqrt(length(soils_PSI$PfracTotAl))*100

ancova(PfracNaOHSRP~PfracNaOHAl*Area, data = soils_PSI)
SRP_NaOH_a <- lm(PfracNaOHSRP~PfracNaOHAl, data = soils_PSI[soils_PSI$Area=="Acidic",])
SRP_NaOH_a_cor <- cor.test(soils_PSI$PfracNaOHSRP[soils_PSI$Area=="Acidic"], soils_PSI$PfracNaOHAl[soils_PSI$Area=="Acidic"])
pvalue_table_pt3[8,2] <- p_value_ext(SRP_NaOH_a, pt = 2, cor = SRP_NaOH_a_cor)

SRP_NaOH_AD <- lm(PfracNaOHSRP~PfracNaOHAl, data = soils_PSI[soils_PSI$treatment=="Acidic Dry",])
SRP_NaOH_AD_cor <- cor.test(soils_PSI$PfracNaOHSRP[soils_PSI$treatment=="Acidic Dry"], soils_PSI$PfracNaOHAl[soils_PSI$treatment=="Acidic Dry"])
pvalue_table_pt3[8,3] <- p_value_ext(SRP_NaOH_AD, pt = 2, cor = SRP_NaOH_AD_cor)

SRP_NaOH_AM <- lm(PfracNaOHSRP~PfracNaOHAl, data = soils_PSI[soils_PSI$treatment=="Acidic Mesic",])
SRP_NaOH_AM_cor <- cor.test(soils_PSI$PfracNaOHSRP[soils_PSI$treatment=="Acidic Mesic"], soils_PSI$PfracNaOHAl[soils_PSI$treatment=="Acidic Mesic"])
pvalue_table_pt3[8,4] <- p_value_ext(SRP_NaOH_AM, pt = 2, cor = SRP_NaOH_AM_cor)

SRP_NaOH_AH <- lm(PfracNaOHSRP~PfracNaOHAl, data = soils_PSI[soils_PSI$treatment=="Acidic Hydric",])
SRP_NaOH_AH_cor <- cor.test(soils_PSI$PfracNaOHSRP[soils_PSI$treatment=="Acidic Hydric"], soils_PSI$PfracNaOHAl[soils_PSI$treatment=="Acidic Hydric"])
pvalue_table_pt3[8,5] <- p_value_ext(SRP_NaOH_AH, pt = 2, cor = SRP_NaOH_AH_cor)

cor.test(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Non-acidic")], soils_PSI$PfracNaOHAl[which(soils_PSI$Area=="Non-acidic")])

# Calcium
Ca_area <- lm(sqrt(PfracHClCa)~Area, data = soils_PSI)
shapiro.test(resid(Ca_area))
anova(Ca_area)
pvalue_table[15,2] <- p_value_ext(Ca_area, pt = 1)

Ca_siteA<- lm(PfracHClCa~SiteName, data = soils_PSI[soils_PSI$Area=="Acidic",])
shapiro.test(resid(Ca_siteA))
anova(Ca_siteA)
base::summary(Ca_siteA)
pvalue_table[15,4] <- p_value_ext(Ca_siteA, pt = 1)

range(soils_PSI$PfracHClCa[which(soils_PSI$Area=="Acidic")])
range(soils_PSI$PfracHClCa[which(soils_PSI$Area=="Non-acidic")])

Ca_siteNA<- lm(PfracHClCa~SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(Ca_siteNA))
anova(Ca_siteNA)
base::summary(Ca_siteNA)
pvalue_table[15,5] <- p_value_ext(Ca_siteNA, pt = 1)

Ca_site<- lm(log(PfracHClCa)~SiteName, data = soils_PSI)
shapiro.test(resid(Ca_site))
anova(Ca_site)
base::summary(Ca_site)
pvalue_table[15,3] <- p_value_ext(Ca_site, pt = 1)

Ca_time <- lm(log(PfracHClCa)~SampleEvent, data = soils_PSI)
shapiro.test(resid(Ca_time))
anova(Ca_time)
base::summary(Ca_time)
pvalue_table[15,6] <- p_value_ext(Ca_time, pt = 1)

mean(soils_PSI$PfracHClCa[soils_PSI$treatment=="Non-acidic Mesic"])
sd(soils_PSI$PfracHClCa[soils_PSI$treatment=="Non-acidic Mesic"])/sqrt(length(soils_PSI$PfracHClCa[soils_PSI$treatment=="Non-acidic Mesic"]))
mean(soils_PSI$PfracHClCa[soils_PSI$treatment=="Non-acidic Hydric"])
sd(soils_PSI$PfracHClCa[soils_PSI$treatment=="Non-acidic Hydric"])/sqrt(length(soils_PSI$PfracHClCa[soils_PSI$treatment=="Non-acidic Hydric"]))

mean(soils_PSI$PfracHClCa[soils_PSI$treatment=="Non-acidic Dry"])
sd(soils_PSI$PfracHClCa[soils_PSI$treatment=="Non-acidic Dry"])/sqrt(length(soils_PSI$PfracHClCa[soils_PSI$treatment=="Non-acidic Dry"]))

TotCa_area <- lm(PfracTotCa~Area, data = soils_PSI)
shapiro.test(resid(TotCa_area))
anova(TotCa_area)
pvalue_table[17,2] <- p_value_ext(TotCa_area, pt = 1)

TotCa_site <- lm(PfracTotCa~SiteName, data = soils_PSI)
shapiro.test(resid(TotCa_site))
anova(TotCa_site)
pvalue_table[17,3] <- p_value_ext(TotCa_site, pt = 1)

TotCa_site_a <- lm(PfracTotCa~SiteName, data = soils_PSI[soils_PSI$Area=="Acidic",])
shapiro.test(resid(TotCa_site_a))
anova(TotCa_site_a)
pvalue_table[17,4] <- p_value_ext(TotCa_site_a, pt = 1)

TotCa_site_na <- lm(PfracTotCa~SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(TotCa_site_na))
anova(TotCa_site_na)
pvalue_table[17,5] <- p_value_ext(TotCa_site_na, pt = 1)

TotCa_time <- lm(PfracTotCa~SampleEvent, data = soils_PSI)
shapiro.test(resid(TotCa_time))
anova(TotCa_time)
pvalue_table[17,6] <- p_value_ext(TotCa_time, pt = 1)

mean(soils_PSI$PfracTotCa[soils_PSI$Area=="Non-acidic"])
sd(soils_PSI$PfracTotCa[soils_PSI$Area=="Non-acidic"])/sqrt(length(soils_PSI$PfracTotCa[soils_PSI$Area=="Non-acidic"]))

mean(soils_PSI$PfracTotCa[soils_PSI$Area=="Acidic"])
sd(soils_PSI$PfracTotCa[soils_PSI$Area=="Acidic"])/sqrt(length(soils_PSI$PfracTotCa[soils_PSI$Area=="Acidic"]))

mean(soils_PSI$PfracHClCa/soils_PSI$PfracTotCa)*100
sd(soils_PSI$PfracHClCa/soils_PSI$PfracTotCa)/sqrt(length(soils_PSI$PfracTotCa))*100

Ca_HCl_area <- lm(lnHClP~PfracHClCa*Area, data = soils_PSI)
Ca_HCl_area_cor <- 
hist(resid(Ca_HCl_area))
anova(Ca_HCl_area)

ancova(lnHClP~PfracHClCa*Area, data = soils_PSI)

Ca_HCl_a <- lm(lnHClP~PfracHClCa*SiteName, data = soils_PSI[which(soils_PSI$Area=="Acidic"),])
shapiro.test(resid(Ca_HCl_E))
anova(Ca_HCl_E)
ancova(lnHClP~PfracHClCa*SiteName, data = soils_PSI[which(soils_PSI$Area=="Acidic"),])

Ca_HCl_na <- lm(lnHClP~PfracHClCa, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(Ca_HCl_na))
anova(Ca_HCl_na)
base::summary(Ca_HCl_na)
Ca_HCl_na_cor <- cor.test(soils_PSI$lnHClP[soils_PSI$Area=="Non-acidic"],soils_PSI$PfracHClCa[soils_PSI$Area=="Non-acidic"])
pvalue_table_pt2[9,2] <- p_value_ext(Ca_HCl_na, pt = 2, Ca_HCl_na_cor)
pvalue_table_pt3[9,2] <- p_value_ext(Ca_HCl_na, pt = 2, Ca_HCl_na_cor)


Ca_HCl_na_site <- lm(lnHClP~PfracHClCa*SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(Ca_HCl_na_site))
anova(Ca_HCl_na_site)
base::summary(Ca_HCl_na_site)
pvalue_table_pt2[9,3] <- p_value_ext(Ca_HCl_na_site, pt = 1)

P_HCLCa_NAD_cor <-cor.test(soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Dry"),]$PfracHClP, soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Dry"),]$PfracHClCa, method = "pearson")
# cor = -0.485, p-value = 0.223
P_HClCa_NAD <- lm(PfracHClP~PfracHClCa, data = soils_PSI[soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Dry",])
pvalue_table_pt3[9,3] <- p_value_ext(P_HClCa_NAD, pt = 2, P_HCLCa_NAD_cor)

P_HClCa_NAM_cor <-cor.test(soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Mesic"),]$PfracHClP, soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Mesic"),]$PfracHClCa, method = "pearson")
# corr = 0.442, p-value = 0.232
P_HClCa_NAM <- lm(PfracHClP~PfracHClCa, data = soils_PSI[soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Mesic",])
pvalue_table_pt3[9,4] <- p_value_ext(P_HClCa_NAM, pt = 2, P_HClCa_NAM_cor)


P_HClCa_NAH_cor <- cor.test(soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Hydric"),]$PfracHClP, soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Hydric"),]$PfracHClCa, method = "pearson")
# cor = 0.527, p-value = 0.144
P_HClCa_NAH <- lm(PfracHClP~PfracHClCa, data = soils_PSI[soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Hydric",])
pvalue_table_pt3[9,5] <- p_value_ext(P_HClCa_NAH, pt =2,P_HClCa_NAH_cor)


### PSI analysis ####
mean(soils_PSI$PfracPSI[which(soils_PSI$SiteName=="Dry")])
sd(soils_PSI$PfracPSI[which(soils_PSI$SiteName=="Dry")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$SiteName=="Dry")]))
mean(soils_PSI$PfracPSI[which(soils_PSI$SiteName=="Mesic")])
sd(soils_PSI$PfracPSI[which(soils_PSI$SiteName=="Mesic")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$SiteName=="Mesic")]))
mean(soils_PSI$PfracPSI[which(soils_PSI$SiteName=="Hydric")])
sd(soils_PSI$PfracPSI[which(soils_PSI$SiteName=="Hydric")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$SiteName=="Hydric")]))

PSI_site <- lm(lnPSI~SiteName, data = soils_PSI)
shapiro.test(resid(PSI_site))
anova(PSI_site)
pvalue_table[18,3] <- p_value_ext(PSI_site, pt =1)

mean(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])
sd(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")]))
mean(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")])
sd(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")]))

PSI_area <- lm(lnPSI~Area, data = soils_PSI)
shapiro.test(resid(PSI_area))
anova(PSI_area)
pvalue_table[18,2] <- p_value_ext(PSI_area, pt = 1)

PSI_site_a <- lm(lnPSI~SiteName, data = soils_PSI[soils_PSI$Area=="Acidic",])
shapiro.test(resid(PSI_site_a))
anova(PSI_site_a)
pvalue_table[18,4] <- p_value_ext(PSI_site_a, pt =1)

PSI_site_na <- lm(lnPSI~SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])
shapiro.test(resid(PSI_site_na))
anova(PSI_site_na)
pvalue_table[18,5] <- p_value_ext(PSI_site_na, pt =1)

PSI_time <- lm(lnPSI~SampleEvent, data = soils_PSI)
shapiro.test(resid(PSI_time))
anova(PSI_time)
pvalue_table[18,6] <- p_value_ext(PSI_time, pt =1)

soils_PSI[which(soils_PSI$PfracPSI==max(soils_PSI$PfracPSI)),1:7]
soils_PSI[which(soils_PSI$PfracPSI==max(soils_PSI$PfracPSI[which(soils_PSI$Area == "Non-acidic")])),1:7]
mean(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic"&soils_PSI$SiteName=="Hydric")])
sd(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic"&soils_PSI$SiteName=="Hydric")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic"&soils_PSI$SiteName=="Hydric")]))

soils_PSI[which(soils_PSI$PfracPSI==min(soils_PSI$PfracPSI)),1:7]
soils_PSI[which(soils_PSI$PfracPSI==min(soils_PSI$PfracPSI[which(soils_PSI$Area == "Non-acidic")])),1:7]
mean(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic"&soils_PSI$SiteName=="Dry")])
sd(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic"&soils_PSI$SiteName=="Dry")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic"&soils_PSI$SiteName=="Dry")]))


H2O_PSI_cor <-cor.test(soils_PSI$PfracPSI, soils_PSI$PfracH2OSRP, method = "spearman")
H2O_PSI <- lm(PfracH2OSRP~PfracPSI, data = soils_PSI)
pvalue_table_pt2[10,2] <- p_value_ext(H2O_PSI, pt = 3,H2O_PSI_cor)

H2O_PSI_DH_site <- lm(PfracH2OSRP~PfracPSI*SiteName, data = soils_PSI)
pvalue_table_pt2[10,3] <- p_value_ext(H2O_PSI_DH_site, pt = 1)

write.csv(pvalue_table, "~/Pfrac Manuscript/anova_table.csv")
write.csv(pvalue_table_pt2, "~/Pfrac Manuscript/ancova_table_pt1.csv")
write.csv(pvalue_table_pt3, "~/Pfrac Manuscript/ancova_table_pt2.csv")

### PSI model selection ####
parameters_na <- soils_PSI %>% filter(Area=="Non-acidic") %>% dplyr::select(LOI,moisture,pH, C.,N., CN_ratio,PfracBDFe,PfracHClCa,PfracNaOHAl,FefracDHFe,sq_Fe,,OrgFe,SiteName) %>%  mutate(SiteName=recode(SiteName, "Dry"=1,"Mesic" = 2,"Hydric"=3))
parm_names <- c("intercept",names(parameters_na))
PSI.pred_na <- regsubsets(x = parameters_na, y = soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")], nbest=4, nvmax = 13,method = c("exhaustive"), df = nrow(parameters_na)-1, weights=rep(1,nrow(parameters_na)))
par(mfrow=c(1,2))
plot(x = PSI.pred_na,labels=parm_names,main = NULL,scale="adjr2",col=gray(seq(0, 0.9, length = 50)))
plot(x = PSI.pred_na,labels=parm_names,main = NULL,scale="bic",col=gray(seq(0, 0.9, length = 50)))

PSI_na_global <-lm(PfracPSI~PfracHClCa+moisture+pH, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
par(mfrow=c(2,2))
plot(PSI_na_global)

PSI.mod_na <- list()
PSI.mod_na[[1]] <- glm(PfracPSI~pH+PfracHClCa+moisture,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[2]] <- glm(PfracPSI~PfracHClCa,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[3]] <- glm(PfracPSI~moisture,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[4]] <- glm(PfracPSI~pH,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[5]] <- glm(PfracPSI~PfracHClCa+moisture,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[6]] <- glm(PfracPSI~PfracHClCa+pH,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[7]] <- glm(PfracPSI~moisture+pH,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[8]] <- glm(PfracPSI~1,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
Modnames <- c("global","Ca","moisture","pH","moisture+Ca","Ca+pH","mositure+pH","intercept")
PSI_na_table <-aictab(cand.set = PSI.mod_na,modnames = Modnames)
confset(cand.set = PSI.mod_na,modnames = Modnames)

modavg(cand.set = PSI.mod_na,parm = "PfracHClCa", modnames = Modnames)
modavg(cand.set = PSI.mod_na,parm = "pH", modnames = Modnames)
modavg(cand.set = PSI.mod_na,parm = "moisture", modnames = Modnames)


PSI_na <- lm(PfracPSI~PfracHClCa+pH+moisture, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
shapiro.test(resid(PSI_na))
anova(PSI_na)
summary(PSI_na)

cor.test(soils_PSI$PfracHClCa[which(soils_PSI$Area=="Non-acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")])
cor.test(soils_PSI$pH[which(soils_PSI$Area=="Non-acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")])
cor.test(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic"&soils_PSI$PfracHClCa <5)], soils_PSI$moisture[which(soils_PSI$Area=="Non-acidic"&soils_PSI$PfracHClCa<5)])
cor.test(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic"&soils_PSI$PfracHClCa <5)], soils_PSI$PfracHClCa[which(soils_PSI$Area=="Non-acidic"&soils_PSI$PfracHClCa<5)])


parameters_a <- soils_PSI %>% filter(Area=="Acidic") %>% dplyr::select(LOI,moisture,pH, C.,N., CN_ratio,PfracBDFe,PfracHClCa,PfracNaOHAl,FefracDHFe,sq_Fe,OrgFe,SiteName) %>% mutate(SiteName=recode(SiteName, "Dry"=1,"Mesic" = 2,"Hydric"=3))
parm_names <- c("intercept",names(parameters_a))
PSI.pred_a <- regsubsets(x = parameters_a, y = soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")], nbest=4, nvmax = 13,method = c("exhaustive"), df = nrow(parameters_a)-1, weights=rep(1,nrow(parameters_a)))
par(mfrow=c(1,2))
plot(x = PSI.pred_a,labels=parm_names,main = NULL,scale="adjr2",col=gray(seq(0, 0.9, length = 50)))
plot(x = PSI.pred_a,labels=parm_names,main = NULL,scale="bic",col=gray(seq(0, 0.9, length = 50)))


PSI_a_global <-lm(PfracPSI~SiteName+PfracNaOHAl+PfracBDFe+moisture, data = soils_PSI[which(soils_PSI$Area=="Acidic"),])
par(mfrow=c(2,2))
plot(PSI_a_global)

PSI.mod_a <- list()
PSI.mod_a[[1]] <- glm(PfracPSI~SiteName+PfracNaOHAl+PfracBDFe+moisture,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[2]] <- glm(PfracPSI~PfracNaOHAl,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[3]] <- glm(PfracPSI~PfracBDFe,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[4]] <- glm(PfracPSI~moisture,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[5]] <- glm(PfracPSI~SiteName,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[6]] <- glm(PfracPSI~PfracNaOHAl+PfracBDFe,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[7]] <- glm(PfracPSI~PfracNaOHAl+PfracBDFe+moisture,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[8]] <- glm(PfracPSI~PfracNaOHAl+PfracBDFe+SiteName,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[9]] <- glm(PfracPSI~PfracBDFe+moisture,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[10]] <- glm(PfracPSI~PfracNaOHAl+moisture,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[11]] <- glm(PfracPSI~PfracBDFe+SiteName,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[12]] <- glm(PfracPSI~PfracNaOHAl+SiteName,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[13]] <- glm(PfracPSI~SiteName+moisture,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
PSI.mod_a[[14]] <- glm(PfracPSI~1,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
Modnames <- c("global","Al","Fe","moisture","SiteName","Al+Fe","Al+Fe+moisture","Al+Fe+SiteName","Fe+moisture","Al+moisture","Fe+SiteName","Al+SiteName","moisture+SiteName","intercept")
PSI_a_table <-aictab(cand.set = PSI.mod_a,modnames = Modnames)
confset(cand.set = PSI.mod_a,modnames = Modnames)
evidence(PSI_a_table, model.high = "Al+Fe+moisture",model.low = "global")
evidence(PSI_a_table, model.high = "Al+Fe",model.low = "global")

modavg(cand.set = PSI.mod_a,parm = "PfracBDFe", modnames = Modnames)
modavg(cand.set = PSI.mod_a,parm = "moisture", modnames = Modnames)
modavg(cand.set = PSI.mod_a,parm = "PfracNaOHAl", modnames = Modnames)

PSI_a <- lm(PfracPSI~PfracBDFe+PfracNaOHAl+moisture, data = soils_PSI[which(soils_PSI$Area=="Acidic"),])
shapiro.test(resid(PSI_a))
anova(PSI_a)
summary(PSI_a)

cor.test(soils_PSI$moisture[which(soils_PSI$Area=="Acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])
cor.test(soils_PSI$PfracBDFe[which(soils_PSI$Area=="Acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])
cor.test(soils_PSI$PfracNaOHAl[which(soils_PSI$Area=="Acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])
cor.test(soils_PSI$sq_Fe[which(soils_PSI$Area=="Acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])
cor.test(soils_PSI$FefracDHFe[which(soils_PSI$Area=="Acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])

Al_moist_PSI <- lm(PfracPSI~PfracNaOHAl*moisture, data = soils_PSI[soils_PSI$Area == "Acidic",])
summary(Al_moist_PSI)
anova(Al_moist_PSI)


HClP.pred_na <- regsubsets(x = parameters_na, y = soils_PSI$PfracHClP[which(soils_PSI$Area=="Non-acidic")], nbest=5, nvmax = 13,method = c("exhaustive"), df = nrow(parameters_na)-1, weights=rep(1,nrow(parameters_na)))
par(mfrow=c(1,2))
plot(x = HClP.pred_na,labels=parm_names,main = NULL,scale="adjr2",col=gray(seq(0, 0.9, length = 50)))
plot(x = HClP.pred_na,labels=parm_names,main = NULL,scale="bic",col=gray(seq(0, 0.9, length = 50)))

HClP_na_global <-lm(PfracHClP~PfracHClCa+PfracBDFe+C.+SiteName, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
par(mfrow=c(2,2))
plot(HClP_na_global)

HClP.mod_na <- list()
HClP.mod_na[[1]] <- glm(PfracHClP~PfracHClCa+PfracBDFe+C.+SiteName,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[2]] <- glm(PfracHClP~SiteName,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[3]] <- glm(PfracHClP~PfracHClCa,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[4]] <- glm(PfracHClP~PfracBDFe,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[5]] <- glm(PfracHClP~C.,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[6]] <- glm(PfracHClP~PfracHClCa+PfracBDFe,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[7]] <- glm(PfracHClP~PfracHClCa+C.,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[8]] <- glm(PfracHClP~PfracHClCa+SiteName,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[9]] <- glm(PfracHClP~PfracBDFe+SiteName,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[10]] <- glm(PfracHClP~PfracBDFe+C.,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[11]] <- glm(PfracHClP~PfracHClCa+C.+SiteName,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[12]] <- glm(PfracHClP~PfracBDFe+C.+SiteName,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
HClP.mod_na[[13]] <- glm(PfracHClP~SiteName+C.,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
HClP.mod_na[[14]] <- glm(PfracHClP~1,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Acidic",])
Modnames <- c("global","Site","Ca","Fe","C%","Ca+Fe","Ca+C%","Ca+Site","Fe+Site","Fe+C%","Ca+C%+Site","Fe+C%+Site","C%+Site","intercept")
HClP_na_table <-aictab(cand.set = HClP.mod_na,modnames = Modnames)
confset(cand.set = HClP.mod_na,modnames = Modnames)
evidence(HClP_na_table, model.high = "C%+Site",model.low = "global")
evidence(HClP_na_table, model.high = "Ca+Fe",model.low = "global")

modavg(cand.set = HClP.mod_na,parm = "PfracHClCa", modnames = Modnames)
modavg(cand.set = HClP.mod_na,parm = "C.", modnames = Modnames)
modavg(cand.set = HClP.mod_na,parm = "SiteName", modnames = Modnames)
