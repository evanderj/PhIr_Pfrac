setwd("~/PhIr/PhIr 2021/For github")
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
shapiro.test(soils_PSI$lnH2OSRP) # still not normal
soils_PSI$lnHClP <- log(soils_PSI$PfracHClP)
shapiro.test(soils_PSI$lnHClP) # normal
soils_PSI$lnPSI <- log(soils_PSI$PfracPSI)
shapiro.test(soils_PSI$lnPSI)
soils_PSI$H2Oinv <- 1/soils_PSI$PfracH2OSRP
soils_PSI$sq_PP <- sqrt(soils_PSI$FefracPPFe)


### summary tables ####
table_stacked <- soils_PSI %>% dplyr::select(treatment,SampleEvent,moisture,depth,PfracPSI,pH,LOI,C.,N.,CN_ratio,PfracBDFe,PfracBDSRP,PfracH2OSRP,PfracHClCa,PfracHClP,PfracNaOHAl,PfracNaOHSRP,PfracNRP,PfracRP,PfracTotAl,PfracTotCa,PfracTotFe,PfracTotP,ResP,FefracDHFe,FefracHHFe,FefracPPFe,ResFe) %>% gather("Fraction", "Value", 3:28)

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
  AV[,2:27] <- round(AV[,2:27],1)
  SE[,2:27] <- round(SE[,2:27],1)
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
  stats_summary$'FefracPPFe' <- paste(AV$'FefracPPFe',"±",SE$'FefracPPFe')
  stats_summary$'PfracBDFe' <- paste(AV$'PfracBDFe',"±",SE$'PfracBDFe')
  stats_summary$'FefracHHFe' <- paste(AV$'FefracHHFe',"±",SE$'FefracHHFe')
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

iron_table <- treatment_summ %>% dplyr::select(treatment, PfracBDFe,FefracPPFe,FefracHHFe,FefracDHFe,ResFe,PfracTotFe) %>% rename("Iron oxides"=PfracBDFe,"Iron~soil organic matter"=FefracPPFe,"Poorly crystalline iron"=FefracHHFe,"Crystalline iron"=FefracDHFe,"Residual iron"=ResFe, "Total Iron"=PfracTotFe)
write.csv(iron_table,"formatted spreadsheets/Iron_table.csv",row.names = F)

minerals_table <- treatment_summ %>% dplyr::select(treatment, PfracNaOHAl,PfracTotAl,PfracHClCa,PfracTotCa,PfracPSI) %>% rename("Aluminum oxides"=PfracNaOHAl, "Total Al"=PfracTotAl,"Calcium carbonate"=PfracHClCa,"Total Ca"=PfracTotCa, "Phosphate sorption index"=PfracPSI)
write.csv(minerals_table,"formatted spreadsheets/minerals_table.csv",row.names = F)


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
shapiro.test(resid(TP_area))
base::summary(TP_area)

# Total P in acidic soils
TP_a_site <- lm(PfracTotP~SiteName, data = soils_PSI[which(soils_PSI$Area == "Acidic"),])
hist(resid(TP_a_site))
anova(TP_a_site)

mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Mesic")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Mesic")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Mesic")]))
mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Dry")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Dry")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Dry")]))
mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Hydric")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Hydric")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Acidic Hydric")]))

# Total P in non-acidic soils
TP_na_site <- lm(PfracTotP~SiteName, data = soils_PSI[which(soils_PSI$Area == "Non-acidic"),])
hist(resid(TP_na_site))
anova(TP_na_site)
base::summary(TP_na_site)

mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Hydric")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Hydric")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Hydric")]))
mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Dry")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Dry")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Dry")]))
mean(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Mesic")], na.rm = TRUE)
sd(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Mesic")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracTotP[which(soils_PSI$treatment=="Non-acidic Mesic")]))

# rP and total P
mean(soils_PSI$PfracRP/soils_PSI$PfracTotP, na.rm = TRUE)*100

# nrP and total P
mean(soils_PSI$PfracNRP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Acidic")], na.rm = TRUE)*100
sd(soils_PSI$PfracNRP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Acidic")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracNRP[which(soils_PSI$Area=="Acidic")]))*100
mean(soils_PSI$PfracNRP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Non-acidic")], na.rm = TRUE)*100
sd(soils_PSI$PfracNRP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Non-acidic")], na.rm = TRUE)/sqrt(length(soils_PSI$PfracNRP[which(soils_PSI$Area=="Non-acidic")]))*100

NRP_check <- soils_PSI %>% dplyr::select(SampleEvent, Area, SiteName, PlotNumber,PfracNRP, PfracTotP)
NRP_check$Perc <- NRP_check$PfracNRP/NRP_check$PfracTotP*100
NRP_perc_area <- lm(Perc~Area, data = NRP_check)
shapiro.test(resid(NRP_perc_area))
anova(NRP_perc_area)

# Res P
mean(soils_PSI$ResP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Acidic")], na.rm = TRUE)*100
sd(soils_PSI$ResP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Acidic")], na.rm = TRUE)/sqrt(length(soils_PSI$ResP[which(soils_PSI$Area=="Acidic")]))*100

mean(soils_PSI$ResP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Non-acidic")], na.rm = TRUE)*100
sd(soils_PSI$ResP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotP[which(soils_PSI$Area=="Non-acidic")], na.rm = TRUE)/sqrt(length(soils_PSI$ResP[which(soils_PSI$Area=="Non-acidic")]))*100


resP_check <- soils_PSI %>% dplyr::select(SampleEvent, Area, SiteName, PlotNumber,ResP, PfracTotP)
resP_check$Perc <- resP_check$ResP/resP_check$PfracTotP*100
resP_perc_area <- lm(Perc~Area, data = resP_check)
shapiro.test(resid(resP_perc_area))
anova(resP_perc_area)


### Reactive Phosphorous fractionation ####
mean(soils_PSI$PfracH2OSRP/soils_PSI$PfracRP)*100
(sd(soils_PSI$PfracH2OSRP/soils_PSI$PfracRP))/sqrt(length(soils_PSI$PfracH2OSRP))*100

## H2OSRP variance ##
H2OSRP_site <- lm(lnH2OSRP~SiteName, data = soils_PSI)
anova(H2OSRP_site)
shapiro.test(resid(H2OSRP_site))
base::summary(H2OSRP_site)

mean(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Mesic")])
sd(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Mesic")])/sqrt(length(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Mesic")]))
mean(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Dry")])
sd(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Dry")])/sqrt(length(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Dry")]))
mean(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Hydric")])
sd(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Hydric")])/sqrt(length(soils_PSI$PfracH2OSRP[which(soils_PSI$SiteName=="Hydric")]))

H2OSRP_area <- lm(lnH2OSRP~Area, data = soils_PSI)
anova(H2OSRP_area)
hist(resid(H2OSRP_area))


## BDSRP variance #
mean(soils_PSI$PfracBDSRP/soils_PSI$PfracRP)*100
(sd(soils_PSI$PfracBDSRP/soils_PSI$PfracRP))/sqrt(length(soils_PSI$PfracBDSRP))*100

mean(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Mesic")])
sd(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Mesic")])/sqrt(length(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Mesic")]))
mean(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Dry")])
sd(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Dry")])/sqrt(length(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Dry")]))
mean(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Hydric")])
sd(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Hydric")])/sqrt(length(soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Hydric")]))

BDrp_site <- lm(PfracBDSRP~SiteName, data = soils_PSI)
shapiro.test(resid(BDrp_site))
anova(BDrp_site)
base::summary(BDrp_site)

BDrp_area <- lm(PfracBDSRP~Area, data = soils_PSI)
anova(BDrp_area)
shapiro.test(resid(BDrp_area))
base::summary(BDrp_area)


## NaOHSRP variance by area and site##
AlSRP_area <- lm(PfracNaOHSRP~Area, data = soils_PSI)
anova(AlSRP_area)
shapiro.test(resid(AlSRP_area))

AlSRP_site <- lm(PfracNaOHSRP~treatment, data = soils_PSI)
anova(AlSRP_site)
shapiro.test(resid(AlSRP_site))
base::summary(AlSRP_site)

mean(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Acidic")])*100
(sd(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Acidic")]))/sqrt(length(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Acidic")]))*100
mean(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Non-acidic")])*100
(sd(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Non-acidic")]))/sqrt(length(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Non-acidic")]))*100


## HClP variance by area and site##
CaP_area <- lm(lnHClP~Area, data = soils_PSI)
anova(CaP_area)
shapiro.test(resid(CaP_area))

CaP_site <- lm(lnHClP~treatment, data = soils_PSI)
anova(CaP_site)
shapiro.test(resid(CaP_site))
base::summary(CaP_site)

mean(soils_PSI$PfracHClP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Acidic")])*100
(sd(soils_PSI$PfracHClP[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Acidic")]))/sqrt(length(soils_PSI$PfracHClP[which(soils_PSI$Area=="Acidic")]))*100
mean(soils_PSI$PfracHClP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Non-acidic")])*100
(sd(soils_PSI$PfracHClP[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracRP[which(soils_PSI$Area=="Non-acidic")]))/sqrt(length(soils_PSI$PfracHClP[which(soils_PSI$Area=="Non-acidic")]))*100


### Iron fractionation ####
Fe_area <- lm(sqrt(PfracTotFe)~Area, data= soils_PSI)
base::summary(Fe_area)
anova(Fe_area)
hist(resid(Fe_area))

Fe_a_site <- lm(sqrt(PfracTotFe)~SiteName, data= soils_PSI[which(soils_PSI$Area=="Acidic"),])
base::summary(Fe_a_site)
anova(Fe_a_site)
shapiro.test(resid(Fe_a_site))

Fe_na_site <- lm(sqrt(PfracTotFe)~SiteName, data= soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
base::summary(Fe_na_site)
anova(Fe_na_site)
shapiro.test(resid(Fe_na_site))

mean(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Mesic")])
sd(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Mesic")])/sqrt(length(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Mesic")]))

mean(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Mesic")])
sd(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Mesic")])/sqrt(length(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Mesic")]))

mean(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Hydric")])
sd(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Hydric")])/sqrt(length(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Acidic Hydric")]))

mean(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Dry")])
sd(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Dry")])/sqrt(length(soils_PSI$PfracTotFe[which(soils_PSI$treatment=="Non-acidic Dry")]))

# iron~soil organic matter
mean(soils_PSI$FefracPPFe[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Acidic")])*100
(sd(soils_PSI$FefracPPFe[which(soils_PSI$Area=="Acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Acidic")])/sqrt(length(soils_PSI$FefracPPFe[which(soils_PSI$Area=="Acidic")])))*100
mean(soils_PSI$FefracPPFe[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Non-acidic")])*100
(sd(soils_PSI$FefracPPFe[which(soils_PSI$Area=="Non-acidic")]/soils_PSI$PfracTotFe[which(soils_PSI$Area=="Non-acidic")])/sqrt(length(soils_PSI$FefracPPFe[which(soils_PSI$Area=="Non-acidic")])))*100

mean(soils_PSI$FefracPPFe[which(soils_PSI$treatment=="Acidic Hydric")])
sd(soils_PSI$FefracPPFe[which(soils_PSI$treatment=="Acidic Hydric")])/sqrt(length(soils_PSI$FefracPPFe[which(soils_PSI$treatment=="Acidic Hydric")]))
mean(soils_PSI$FefracPPFe[which(soils_PSI$treatment=="Acidic Dry")])
sd(soils_PSI$FefracPPFe[which(soils_PSI$treatment=="Acidic Dry")])/sqrt(length(soils_PSI$FefracPPFe[which(soils_PSI$treatment=="Acidic Dry")]))
mean(soils_PSI$FefracPPFe[which(soils_PSI$treatment=="Acidic Mesic")])
sd(soils_PSI$FefracPPFe[which(soils_PSI$treatment=="Acidic Mesic")])/sqrt(length(soils_PSI$FefracPPFe[which(soils_PSI$treatment=="Acidic Mesic")]))


PP_Fe_a_site <- lm(sqrt(FefracPPFe)~SiteName, data = soils_PSI[which(soils_PSI$Area == "Acidic"),])
base::summary(PP_Fe_a_site)
shapiro.test(resid(PP_Fe_a_site))

PP_Fe_na_site <- lm(sqrt(FefracPPFe)~SiteName, data = soils_PSI[which(soils_PSI$Area == "Non-acidic"),])
base::summary(PP_Fe_na_site)
shapiro.test(resid(PP_Fe_na_site))


# poorly crystalline iron
HH_Fe_area <- lm(sqrt(FefracHHFe)~Area, data = soils_PSI)
shapiro.test(resid(HH_Fe_area))
base::summary(HH_Fe_area)

mean(soils_PSI$FefracHHFe[soils_PSI$Area=="Non-acidic"])
sd(soils_PSI$FefracHHFe[soils_PSI$Area=="Non-acidic"])/sqrt(length(soils_PSI$FefracHHFe[soils_PSI$Area=="Non-acidic"]))
mean(soils_PSI$FefracHHFe[soils_PSI$Area=="Acidic"])
sd(soils_PSI$FefracHHFe[soils_PSI$Area=="Acidic"])/sqrt(length(soils_PSI$FefracHHFe[soils_PSI$Area=="Acidic"]))


HH_Fe_a_site <- lm(sqrt(FefracHHFe)~SiteName, data = soils_PSI[which(soils_PSI$Area=="Acidic"),])
shapiro.test(resid(HH_Fe_a_site))
base::summary(HH_Fe_a_site)

HH_Fe_na_site <- lm(sqrt(FefracHHFe)~SiteName, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
shapiro.test(resid(HH_Fe_na_site))
base::summary(HH_Fe_na_site)

# poorly crystalline iron and iron ~soil organic matter
cor.test(soils_PSI$FefracPPFe[soils_PSI$SiteName=="Hydric"],soils_PSI$FefracHHFe[soils_PSI$SiteName=="Hydric"])

# iron oxides and iron~soil organic matter
cor.test(log(soils_PSI$PfracBDFe), log(soils_PSI$FefracPPFe), method = "pearson")
plot(log(soils_PSI$PfracBDFe), log(soils_PSI$FefracPPFe))
BD_PP <- lm(log(PfracBDFe)~ log(FefracPPFe), data = soils_PSI)
shapiro.test(resid(BD_PP))

# iron oxides and crystalline iron
cor.test(log(soils_PSI$PfracBDFe), log(soils_PSI$FefracDHFe), method = "pearson")
plot(log(soils_PSI$PfracBDFe), log(soils_PSI$FefracDHFe))
BD_DH <- lm(log(PfracBDFe)~ log(FefracDHFe), data = soils_PSI)
shapiro.test(resid(BD_DH))

# iron oxides and poorly crystalline iron
cor.test(log(soils_PSI$PfracBDFe), soils_PSI$FefracHHFe, method = "pearson")
plot(log(soils_PSI$PfracBDFe), soils_PSI$FefracHHFe)
BD_HH <- lm(log(PfracBDFe)~ FefracHHFe, data = soils_PSI)
shapiro.test(resid(BD_HH))


### reactive phosphorus pools ####
## BDSRP covariance with H2OSRP ##
SRP_H2O <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI)
anova(SRP_H2O)
base::summary(SRP_H2O)
shapiro.test(resid(SRP_H2O))

ancova(lnH2OSRP~PfracBDSRP*SiteName, data = soils_PSI)
SRP_H2O_site <- lm(lnH2OSRP~PfracBDSRP*SiteName, data = soils_PSI)
base::summary(SRP_H2O_site)

SRP_H2O_D <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI[which(soils_PSI$SiteName=="Dry"),])
base::summary(SRP_H2O_D)
coef(SRP_H2O_D)
cor.test(soils_PSI$lnH2OSRP[which(soils_PSI$SiteName=="Dry")], soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Dry")])

SRP_H2O_M <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI[which(soils_PSI$SiteName=="Mesic"),])
base::summary(SRP_H2O_M)
coef(SRP_H2O_M)
cor.test(soils_PSI$lnH2OSRP[which(soils_PSI$SiteName=="Mesic")], soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Mesic")])

SRP_H2O_H <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI[which(soils_PSI$SiteName=="Hydric"),])
base::summary(SRP_H2O_H)
coef(SRP_H2O_H)
cor.test(soils_PSI$lnH2OSRP[which(soils_PSI$SiteName=="Hydric")], soils_PSI$PfracBDSRP[which(soils_PSI$SiteName=="Hydric")])


## SRP covariance between H2O, NaOH, and HCl ##

ancova(lnH2OSRP~PfracNaOHSRP*SiteName, data = soils_PSI)
NaOH_H2O_site <- lm(lnH2OSRP~PfracNaOHSRP, data = soils_PSI[which(soils_PSI$SiteName!="Dry"),])
base::summary(NaOH_H2O_site)

NaOH_H2O_M <- lm(lnH2OSRP~PfracNaOHSRP, data = soils_PSI[which(soils_PSI$SiteName=="Mesic"),])
base::summary(NaOH_H2O_M)
cor.test(soils_PSI$PfracNaOHSRP[which(soils_PSI$SiteName == "Mesic")], soils_PSI$lnH2OSRP[which(soils_PSI$SiteName == "Mesic")], method = "pearson")
NaOH_H2O_H <- lm(lnH2OSRP~PfracNaOHSRP, data = soils_PSI[which(soils_PSI$SiteName=="Hydric"),])
base::summary(NaOH_H2O_H)
cor.test(soils_PSI$PfracNaOHSRP[which(soils_PSI$SiteName == "Hydric")], soils_PSI$lnH2OSRP[which(soils_PSI$SiteName == "Hydric")], method = "pearson")


ancova(lnH2OSRP~PfracHClP*SiteName, data = soils_PSI)
ancova(lnH2OSRP~PfracHClP*Area, data = soils_PSI)
ancova(lnH2OSRP~PfracHClP*treatment, data = soils_PSI)


## BDSRP and NaOHSRP, HClP ##
SRP_BD_NaOH <- lm(PfracBDSRP~PfracNaOHSRP*Area, data = soils_PSI)
anova(SRP_BD_NaOH)
base::summary(SRP_BD_NaOH)
shapiro.test(resid(SRP_BD_NaOH))
ancova(PfracBDSRP~PfracNaOHSRP*Area, data = soils_PSI)

cor.test(soils_PSI$PfracBDSRP[which(soils_PSI$Area == "Non-acidic")], soils_PSI$PfracNaOHSRP[which(soils_PSI$Area == "Non-acidic")], method = "pearson")
cor.test(soils_PSI$PfracBDSRP[which(soils_PSI$Area == "Acidic")], soils_PSI$PfracNaOHSRP[which(soils_PSI$Area == "Acidic")], method = "pearson")

SRP_BD_HCl <- lm(PfracBDSRP~PfracHClP, data = soils_PSI)
anova(SRP_BD_HCl)
shapiro.test(resid(SRP_BD_HCl))

SRP_BD_HCl_area <- lm(PfracBDSRP~PfracHClP*Area, data = soils_PSI)
anova(SRP_BD_HCl_area)
shapiro.test(resid(SRP_BD_HCl_area))

SRP_BD_HCl_site <- lm(PfracBDSRP~PfracHClP*SiteName, data = soils_PSI)
anova(SRP_BD_HCl_site)
shapiro.test(resid(SRP_BD_HCl_site))


### H2OSRP and minerals ####
## H2OSRP and iron oxides##
H2Orp_BDFe <- lm(H2Oinv~PfracBDFe, data = soils_PSI)
hist(resid(H2Orp_BDFe))
anova(H2Orp_BDFe)

ancova(H2Oinv~PfracBDFe*treatment, data = soils_PSI)
cor.test(soils_PSI$PfracH2OSRP, soils_PSI$PfracBDFe)
# r = -0.27, p = 0.05

# H2OSRP and iron~soil organic matter
H2Orp_PPFe <- lm(H2Oinv~sqrt(FefracPPFe), data = soils_PSI)
shapiro.test(resid(H2Orp_PPFe))
anova(H2Orp_PPFe)

cor.test(soils_PSI$PfracH2OSRP, soils_PSI$FefracPPFe)
# r = -0.21, p = 0.13

# H2OSRP and poorly crystalline iron
H2Orp_HHFe <- lm(lnH2OSRP~FefracHHFe, data = soils_PSI)
shapiro.test(resid(H2Orp_HHFe))
anova(H2Orp_HHFe)

cor.test(soils_PSI$PfracH2OSRP, soils_PSI$FefracHHFe, method = "pearson")
# r = 0.24, p = 0.09

# H2OSRP and crystalline iron
H2Orp_DHFe <- lm(H2Oinv~FefracDHFe, data = soils_PSI)
shapiro.test(resid(H2Orp_DHFe))
anova(H2Orp_DHFe)
cor.test(soils_PSI$PfracH2OSRP, soils_PSI$FefracDHFe, method = "pearson")
# r = -0.1, p = 0.46


## H2O~SRP vs Ca and Al ##
H2O_Al <- lm(H2Oinv~PfracNaOHAl, data = soils_PSI)
shapiro.test(resid(H2O_Al))
anova(H2O_Al)
base::summary(H2O_Al)

cor.test(soils_PSI$PfracH2OSRP, soils_PSI$PfracNaOHAl, method = "pearson")

H2O_Ca <- lm(lnH2OSRP~PfracHClCa*Area, data = soils_PSI)
shapiro.test(resid(H2O_Ca))
ancova(PfracH2OSRP~PfracHClCa*Area, data = soils_PSI)
base::summary(H2O_Ca)

ancova(PfracH2OSRP~PfracHClCa*SiteName, data = soils_PSI[soils_PSI$Area=="Non-acidic",])

cor.test(soils_PSI$PfracH2OSRP[which(soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Dry")],soils_PSI$PfracHClCa[which(soils_PSI$Area == "Non-acidic" & soils_PSI$SiteName=="Dry")], method = "pearson")


### BD~SRP and iron ####
## BD~SRP and iron~soil organic matter ##
PPFe_BDSRP_area <- lm(PfracBDSRP~FefracPPFe+Area, data = soils_PSI)
hist(resid(PPFe_BDSRP_area))
anova(PPFe_BDSRP_area)

PPFe_BDSRP_site <- lm(PfracBDSRP~FefracPPFe*SiteName, data = soils_PSI)
hist(resid(PPFe_BDSRP_site))
anova(PPFe_BDSRP_site)
base::summary(PPFe_BDSRP_site)

PPFe_BDSRP_treat <- lm(PfracBDSRP~FefracPPFe+treatment, data = soils_PSI)
hist(resid(PPFe_BDSRP_treat))
anova(PPFe_BDSRP_treat)

ancova(PfracBDSRP~sq_PP*treatment, data = soils_PSI)

cor.test(soils_PSI[which(soils_PSI$SiteName=="Hydric"& soils_PSI$Area =="Acidic"),]$PfracBDSRP, soils_PSI[which(soils_PSI$SiteName=="Hydric"& soils_PSI$Area =="Acidic"),]$FefracPPFe, method = "pearson")
# corr = 0.467, p-value = 0.205
cor.test(soils_PSI[which(soils_PSI$SiteName=="Hydric"& soils_PSI$Area =="Non-acidic"),]$PfracBDSRP, soils_PSI[which(soils_PSI$SiteName=="Hydric"& soils_PSI$Area =="Non-acidic"),]$FefracPPFe, method = "pearson")
# corr = 0.63, p-value = 0.07

cor.test(soils_PSI[which(soils_PSI$SiteName=="Dry"),]$PfracBDSRP, log(soils_PSI[which(soils_PSI$SiteName=="Dry"),]$FefracPPFe), method = "pearson")

ancova(PfracBDSRP~FefracHHFe*treatment, data = soils_PSI)
ancova(PfracBDSRP~FefracDHFe*treatment, data = soils_PSI)

cor.test(soils_PSI[which(soils_PSI$treatment=="Acidic Dry"),]$PfracBDSRP, soils_PSI[which(soils_PSI$treatment=="Acidic Dry"),]$FefracDHFe, method = "pearson")
cor.test(soils_PSI[which(soils_PSI$treatment=="Acidic Mesic"),]$PfracBDSRP, soils_PSI[which(soils_PSI$treatment=="Acidic Mesic"),]$FefracDHFe, method = "pearson")


## Aluminum and Calcium  ##
Al_area <- lm(PfracNaOHAl~Area, soils_PSI)
shapiro.test(resid(Al_area))
anova(Al_area)

range(soils_PSI$PfracNaOHAl[which(soils_PSI$Area=="Acidic")])
range(soils_PSI$PfracNaOHAl[which(soils_PSI$Area=="Non-acidic")])

ancova(PfracNaOHSRP~PfracNaOHAl*Area, data = soils_PSI)
cor.test(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Acidic")], soils_PSI$PfracNaOHAl[which(soils_PSI$Area=="Acidic")])

cor.test(soils_PSI$PfracNaOHSRP[which(soils_PSI$Area=="Non-acidic")], soils_PSI$PfracNaOHAl[which(soils_PSI$Area=="Non-acidic")])


Ca_area <- lm(PfracHClCa~Area, data = soils_PSI)
hist(resid(Ca_area))
anova(Ca_area)

range(soils_PSI$PfracHClCa[which(soils_PSI$Area=="Acidic")])
range(soils_PSI$PfracHClCa[which(soils_PSI$Area=="Non-acidic")])


Ca_HCl_area <- lm(lnHClP~PfracHClCa*Area, data = soils_PSI)
hist(resid(Ca_HCl_area))
anova(Ca_HCl_area)

ancova(lnHClP~PfracHClCa*Area, data = soils_PSI)

Ca_HCl_E <- lm(lnHClP~PfracHClCa*SiteName, data = soils_PSI[which(soils_PSI$Area=="Acidic"),])
shapiro.test(resid(Ca_HCl_E))
anova(Ca_HCl_E)
ancova(lnHClP~PfracHClCa*SiteName, data = soils_PSI[which(soils_PSI$Area=="Acidic"),])

Ca_HCl_W <- lm(PfracHClP~PfracHClCa*SiteName, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
shapiro.test(resid(Ca_HCl_W))
anova(Ca_HCl_W)
summary(Ca_HCl_W)
ancova(PfracHClP~PfracHClCa*SiteName, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])


cor.test(soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Dry"),]$PfracHClP, soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Dry"),]$PfracHClCa, method = "pearson")
# cor = -0.485, p-value = 0.223

cor.test(soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Mesic"),]$PfracHClP, soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Mesic"),]$PfracHClCa, method = "pearson")
# corr = 0.442, p-value = 0.232

cor.test(soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Hydric"),]$PfracHClP, soils_PSI[which(soils_PSI$Area=="Non-acidic" & soils_PSI$SiteName=="Hydric"),]$PfracHClCa, method = "pearson")
# cor = 0.527, p-value = 0.144


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


mean(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])
sd(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")]))
mean(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")])
sd(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")]))

PSI_date <- lm(lnPSI~Area, data = soils_PSI)
shapiro.test(resid(PSI_date))
anova(PSI_date)

soils_PSI[which(soils_PSI$PfracPSI==max(soils_PSI$PfracPSI)),1:7]
soils_PSI[which(soils_PSI$PfracPSI==max(soils_PSI$PfracPSI[which(soils_PSI$Area == "Non-acidic")])),1:7]
mean(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic"&soils_PSI$SiteName=="Hydric")])
sd(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic"&soils_PSI$SiteName=="Hydric")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic"&soils_PSI$SiteName=="Hydric")]))

soils_PSI[which(soils_PSI$PfracPSI==min(soils_PSI$PfracPSI)),1:7]
soils_PSI[which(soils_PSI$PfracPSI==min(soils_PSI$PfracPSI[which(soils_PSI$Area == "Non-acidic")])),1:7]
mean(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic"&soils_PSI$SiteName=="Dry")])
sd(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic"&soils_PSI$SiteName=="Dry")])/sqrt(length(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic"&soils_PSI$SiteName=="Dry")]))


cor.test(soils_PSI$PfracPSI[soils_PSI$SiteName!="Mesic"], soils_PSI$PfracH2OSRP[soils_PSI$SiteName!="Mesic"])

### PSI model selection ####
parameters_na <- soils_PSI %>% filter(Area=="Non-acidic") %>% select(LOI,moisture,pH, C.,N., CN_ratio,PfracBDFe,PfracHClCa,PfracNaOHAl,FefracDHFe,FefracHHFe,FefracPPFe,SiteName) %>%  mutate(SiteName=recode(SiteName, "Dry"=1,"Mesic" = 2,"Hydric"=3))
parm_names <- c("intercept",names(parameters_na))
PSI.pred_na <- regsubsets(x = parameters_na, y = soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")], nbest=4, nvmax = 13,method = c("exhaustive"), df = nrow(parameters_na)-1, weights=rep(1,nrow(parameters_na)))
par(mfrow=c(1,2))
plot(x = PSI.pred_na,labels=parm_names,main = NULL,scale="adjr2",col=gray(seq(0, 0.9, length = 50)))
plot(x = PSI.pred_na,labels=parm_names,main = NULL,scale="bic",col=gray(seq(0, 0.9, length = 50)))

PSI_na_global <-lm(PfracPSI~PfracHClCa+FefracHHFe+pH, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
par(mfrow=c(2,2))
plot(PSI_na_global)

PSI.mod_na <- list()
PSI.mod_na[[1]] <- glm(PfracPSI~pH+PfracHClCa+FefracHHFe,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[2]] <- glm(PfracPSI~PfracHClCa,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[3]] <- glm(PfracPSI~FefracHHFe,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[4]] <- glm(PfracPSI~pH,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[5]] <- glm(PfracPSI~PfracHClCa+FefracHHFe,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[6]] <- glm(PfracPSI~PfracHClCa+pH,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[7]] <- glm(PfracPSI~FefracHHFe+pH,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
PSI.mod_na[[8]] <- glm(PfracPSI~1,family ="gaussian",data=soils_PSI[soils_PSI$Area=="Non-acidic",])
Modnames <- c("global","Ca","HHFe","pH","HHFe+Ca","Ca+pH","HHFe+pH","intercept")
PSI_na_table <-aictab(cand.set = PSI.mod_na,modnames = Modnames)
confset(cand.set = PSI.mod_na,modnames = Modnames)

modavg(cand.set = PSI.mod_na,parm = "PfracHClCa", modnames = Modnames)
modavg(cand.set = PSI.mod_na,parm = "pH", modnames = Modnames)
modavg(cand.set = PSI.mod_na,parm = "FefracHHFe", modnames = Modnames)


PSI_na <- lm(PfracPSI~PfracHClCa+pH+FefracHHFe, data = soils_PSI[which(soils_PSI$Area=="Non-acidic"),])
shapiro.test(resid(PSI_na))
anova(PSI_na)
summary(PSI_na)

cor.test(soils_PSI$PfracHClCa[which(soils_PSI$Area=="Non-acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")])
cor.test(soils_PSI$pH[which(soils_PSI$Area=="Non-acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic")])
cor.test(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic"&soils_PSI$PfracHClCa <5)], soils_PSI$FefracHHFe[which(soils_PSI$Area=="Non-acidic"&soils_PSI$PfracHClCa<5)])
cor.test(soils_PSI$PfracPSI[which(soils_PSI$Area=="Non-acidic"&soils_PSI$PfracHClCa <5)], soils_PSI$PfracHClCa[which(soils_PSI$Area=="Non-acidic"&soils_PSI$PfracHClCa<5)])


parameters_a <- soils_PSI %>% filter(Area=="Acidic") %>% dplyr::select(LOI,moisture,pH, C.,N., CN_ratio,PfracBDFe,PfracHClCa,PfracNaOHAl,FefracDHFe,FefracHHFe,FefracPPFe,SiteName) %>% mutate(SiteName=recode(SiteName, "Dry"=1,"Mesic" = 2,"Hydric"=3))
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
evidence(PSI_a_table, model.high = "PPFe+PPFe+moisture",model.low = "global")

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
cor.test(soils_PSI$sq_PP[which(soils_PSI$Area=="Acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])
cor.test(soils_PSI$FefracHHFe[which(soils_PSI$Area=="Acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])
cor.test(soils_PSI$FefracDHFe[which(soils_PSI$Area=="Acidic")], soils_PSI$PfracPSI[which(soils_PSI$Area=="Acidic")])

Al_moist_PSI <- lm(PfracPSI~PfracNaOHAl*moisture, data = soils_PSI[soils_PSI$Area == "Acidic",])
summary(Al_moist_PSI)
anova(Al_moist_PSI)
