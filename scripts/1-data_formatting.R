setwd("~/PhIr/PhIr 2021/for github")
source("0-packages.R")

### Phosphorous fractionation data ####
Pfrac<-read.csv("raw data/PhIr2021_PfracSummary.csv", na.strings=c("#DIV/0!", "#VALUE!", " ", "NA"))
PfracFieldCodes<-read.csv("raw data/PFractionation_FieldCodes.csv")
colnames(Pfrac)<-PfracFieldCodes$RColumnHeaders
Pfrac<-Pfrac %>% mutate(SiteName=factor(SiteName, levels = c("Dry", "Mesic", "Hydric"))) %>% mutate(Horizon = factor(Horizon, levels = c("Organic", "Mineral"))) %>% mutate(Area = recode(Area, "West"="Non-acidic", "East"="Acidic"))

##Identify columns w/ negative values
##Replace below detection limit values (For now, replace negatives w/ 0)
Pfrac$PfracH2ONRP[Pfrac$PfracH2ONRP<0]<- 0
Pfrac$PfracH2OCa[Pfrac$PfracH2OCa<0]<- 0
Pfrac$PfracBDNRP[Pfrac$PfracBDNRP<0]<-0

## remove erroneous samples 
Pfrac <- Pfrac[-which(Pfrac$Area=="Acidic"&Pfrac$SiteName=="Mesic"&Pfrac$PlotNumber==1& Pfrac$SoilSampleDate=="2-Jul-21"),]
Pfrac <- Pfrac[-which(Pfrac$Area=="Non-acidic"&Pfrac$SiteName=="Dry"&Pfrac$PlotNumber==2& Pfrac$SoilSampleDate=="22-Jul-21"),]
Pfrac <- Pfrac[-which(Pfrac$Area=="Acidic"&Pfrac$SiteName=="Mesic"&Pfrac$PlotNumber==2& Pfrac$SoilSampleDate=="2-Jul-21"),]
Pfrac <- Pfrac[-which(Pfrac$Area=="Acidic"&Pfrac$SiteName=="Mesic"&Pfrac$PlotNumber==2& Pfrac$SoilSampleDate=="19-Jul-21"),]

#calculating NRP, RP, and ResP 
Pfrac<-Pfrac %>%
  mutate(PfracNRP=PfracH2ONRP+PfracBDNRP+PfracHAP+PfracNaOHNRP,
         PfracRP=PfracH2OSRP+PfracBDSRP+PfracNaOHSRP+PfracHClP) %>%
  mutate(ResP = PfracTotP - PfracRP - PfracNRP,
         PfracNaOHNRP = PfracNaOHNRP + PfracHAP)
names(Pfrac)


##calculate averages for each variable from analytical replicates
Pfrac_Stacked<-gather(Pfrac, "Fraction", "Value", 13:42)
Pfrac_Stacked_Summary<-Pfrac_Stacked %>%
  group_by(SoilSampleDate, Area, SiteName, Fraction, Horizon,PlotNumber)%>%
  summarize(mean=mean(Value, na.rm = TRUE),
            n=length(Value),
            sd=sd(Value, na.rm = TRUE),
            se=(sd(Value, na.rm = TRUE)/sqrt(length(Value)) ),
            rsd=(sd(Value, na.rm = TRUE)/mean(Value))*100
  )


Pfrac_Means<-Pfrac_Stacked_Summary %>%
  dplyr::select(SoilSampleDate, Area, SiteName, PlotNumber,Fraction, Horizon, mean) %>%
  spread(Fraction, mean)

Pfrac_se<-Pfrac_Stacked_Summary %>%
  dplyr::select(SoilSampleDate, Area, SiteName, PlotNumber,Fraction, Horizon, se) %>%
  spread(Fraction, se)


# convert samples dates to sample events
Pfrac<- Pfrac %>% mutate(SampleEvent = recode(SoilSampleDate, "2-Jul-21"="Early Season", "5-Jul-21"="Early Season","19-Jul-21"="Mid Season","22-Jul-21"="Mid Season","2-Aug-21"="Late Season","5-Aug-21"="Late Season"))

Pfrac_Means <- Pfrac_Means %>% mutate(SampleEvent = recode(SoilSampleDate, "2-Jul-21"="Early Season", "5-Jul-21"="Early Season","19-Jul-21"="Mid Season","22-Jul-21"="Mid Season","2-Aug-21"="Late Season","5-Aug-21"="Late Season"))

Pfrac_se <- Pfrac_se %>% mutate(SampleEvent = recode(SoilSampleDate, "2-Jul-21"="Early Season", "5-Jul-21"="Early Season","19-Jul-21"="Mid Season","22-Jul-21"="Mid Season","2-Aug-21"="Late Season","5-Aug-21"="Late Season"))


### Stacked bar graphs ####
Pfrac_Means <- Pfrac_Means[-which(Pfrac_Means$Horizon=="Mineral"),]
Pfrac_Means_stacked <- Pfrac_Means %>% gather("Fraction","Value", 6:35)
Pfrac_snap<-Pfrac_Means_stacked %>%
  group_by(Area, SiteName, Fraction, Horizon)%>%
  summarize(mean=mean(Value, na.rm = TRUE),
            se=(sd(Value, na.rm = TRUE)/sqrt(length(Value))))


# adding together average values to stack values and add error bars to the figures
Pfrac_snap_RP<-Pfrac_snap %>% filter(Fraction == "PfracBDSRP"|Fraction =="PfracH2OSRP"|Fraction == "PfracNaOHSRP"|Fraction == "PfracHClP") %>% group_by(Area,SiteName) %>%
  mutate(Accumulation = if_else(Fraction =="PfracH2OSRP",mean[Fraction=="PfracH2OSRP"]+mean[Fraction=="PfracBDSRP"]+mean[Fraction=="PfracNaOHSRP"]+mean[Fraction=="PfracHClP"],ifelse(Fraction =="PfracBDSRP",mean[Fraction=="PfracBDSRP"]+mean[Fraction=="PfracNaOHSRP"]+mean[Fraction=="PfracHClP"],ifelse(Fraction =="PfracNaOHSRP",mean[Fraction=="PfracNaOHSRP"]+mean[Fraction=="PfracHClP"],ifelse(Fraction=="PfracHClP", mean[Fraction=="PfracHClP"],NA)))))


Pfrac_snap_TP <- Pfrac_snap %>% filter(Fraction == "PfracRP"|Fraction =="PfracNRP"|Fraction == "ResP") %>% group_by(Area,SiteName) %>%
  mutate(Accumulation = if_else(Fraction =="PfracRP",mean[Fraction=="PfracRP"]+mean[Fraction=="PfracNRP"]+mean[Fraction=="ResP"],ifelse(Fraction =="PfracNRP",mean[Fraction=="PfracNRP"]+mean[Fraction=="ResP"],ifelse(Fraction =="ResP",mean[Fraction=="ResP"],NA))))

write.csv(Pfrac_snap_RP, "formatted spreadsheets/Pfrac_RP_stacked.csv",row.names = F)
write.csv(Pfrac_snap_TP, "formatted spreadsheets/Pfrac_TP_stacked.csv", row.names = F)



#### Iron fractionation data ###
Fefrac<-read.csv("raw data/PhIr2021_FeFracSummary.csv", na.strings=c("#DIV/0!", "#VALUE!", "NA"))
FefracFieldCodes<-read.csv("raw data/FeFractionation_FieldCodes.csv")
colnames(Fefrac)<-FefracFieldCodes$RColumnHeaders
Fefrac <-Fefrac %>% mutate(SiteName =factor(SiteName, levels = c("Dry", "Mesic", "Hydric"))) %>% mutate(Horizon = factor(Horizon,  levels = c("Organic", "Mineral"))) %>% mutate(Area=recode(Area,"East"="Acidic","West"="Non-acidic"))

##Calculate Residual Iron
Fefrac$ResFe <- Fefrac$TotFe - Fefrac$FefracPPFe-Fefrac$FefracHHFe-Fefrac$FefracDHFe

# remove erroneous samples
Fefrac <- Fefrac[-which(Fefrac$Area=="Acidic"&Fefrac$SiteName=="Mesic"&Fefrac$PlotNumber==1& Fefrac$SoilSampleDate=="2-Jul-21"),]
Fefrac <- Fefrac[-which(Fefrac$Area=="Non-acidic"&Fefrac$SiteName=="Dry"&Fefrac$PlotNumber==2& Fefrac$SoilSampleDate=="22-Jul-21"),]
Fefrac <- Fefrac[-which(Fefrac$Area=="Acidic"&Fefrac$SiteName=="Mesic"&Fefrac$PlotNumber==2& Fefrac$SoilSampleDate=="2-Jul-21"),]
Fefrac <- Fefrac[-which(Fefrac$Area=="Acidic"&Fefrac$SiteName=="Mesic"&Fefrac$PlotNumber==2& Fefrac$SoilSampleDate=="19-Jul-21"),]


##calculate averages for each variable from analytical replicates
Fefrac_Stacked<-gather(Fefrac, "Fraction", "Value", 13:21)
Fefrac_Stacked_Summary<-Fefrac_Stacked %>%
  group_by(SoilSampleDate,Area,SiteName,Fraction,Horizon,PlotNumber)%>%
  summarize(mean=mean(Value, na.rm= TRUE),
            n=length(Value),
            sd=sd(Value, na.rm = TRUE),
            se=(sd(Value, na.rm = TRUE)/sqrt(length(Value)) ),
            rsd=(sd(Value, na.rm = TRUE)/mean(Value, na.rm = TRUE))*100
  )


##Generating table with QA QC parameters
Fefrac_Means<-Fefrac_Stacked_Summary %>%
  dplyr::select(SoilSampleDate,Area,SiteName,PlotNumber,Fraction,Horizon, mean) %>%
  spread(Fraction, mean)


#se
Fefrac_se<-Fefrac_Stacked_Summary %>%
  dplyr::select(SoilSampleDate,Area,SiteName,PlotNumber,Fraction,Horizon, se) %>%
  spread(Fraction, se)

# convert samples dates to sample events
Fefrac<- Fefrac %>% mutate(SampleEvent = recode(SoilSampleDate, "2-Jul-21"="Early Season", "5-Jul-21"="Early Season","19-Jul-21"="Mid Season","22-Jul-21"="Mid Season","2-Aug-21"="Late Season","5-Aug-21"="Late Season"))

Fefrac_Means <- Fefrac_Means %>% mutate(SampleEvent = recode(SoilSampleDate, "2-Jul-21"="Early Season", "5-Jul-21"="Early Season","19-Jul-21"="Mid Season","22-Jul-21"="Mid Season","2-Aug-21"="Late Season","5-Aug-21"="Late Season"))

Fefrac_se <- Fefrac_se %>% mutate(SampleEvent = recode(SoilSampleDate, "2-Jul-21"="Early Season", "5-Jul-21"="Early Season","19-Jul-21"="Mid Season","22-Jul-21"="Mid Season","2-Aug-21"="Late Season","5-Aug-21"="Late Season"))


### Stacked bar graphs ####
Fefrac_Means <- Fefrac_Means[-which(Fefrac_Means$Horizon == "Mineral"),]
Fefrac_Means_stacked<-gather(Fefrac_Means, "Fraction", "Value", 6:14)
Fefrac_snap<-Fefrac_Means_stacked %>%
  group_by(Area,SiteName,Fraction,Horizon)%>%
  summarize(mean=mean(Value, na.rm= TRUE),
            se=(sd(Value, na.rm = TRUE)/sqrt(length(Value))))


# adding together average values to stack values and add error bars to the figures
Fefrac_snap_Fe<-Fefrac_snap %>% filter(Fraction == "FefracPPFe"|Fraction =="FefracHHFe"|Fraction == "FefracDHFe") %>% group_by(Area,SiteName) %>%
  mutate(Accumulation = if_else(Fraction =="FefracPPFe",mean[Fraction=="FefracPPFe"]+mean[Fraction=="FefracHHFe"]+mean[Fraction=="FefracDHFe"],ifelse(Fraction =="FefracHHFe",mean[Fraction=="FefracHHFe"]+mean[Fraction=="FefracDHFe"],ifelse(Fraction =="FefracDHFe",mean[Fraction=="FefracDHFe"],NA))))

write.csv(Fefrac_snap_Fe, "formatted spreadsheets/Fefrac_Fe_stacked.csv", row.names = F)


### soil conditions ####
wwdw <- read.csv("raw data/PhIr2021_moisture.csv")
LOI <- read.csv("raw data/PhIr2021_LOI.csv")
pH <- read.csv("raw data/soil_pH.csv")
CN <- read.csv("raw data/PhIr2021_CN.csv")

wwdw<-wwdw %>% dplyr::select(Sample.Name,collection.date,Area,Site,Plot,Rep,Horizon,dw.ww,Moisture)
LOI <- LOI %>% dplyr::select(Sample.Name,collection.date,Area,Site,Plot,Rep,Horizon,LOI)
pH <- pH %>% dplyr::select(collection.date,Area,Site,Plot,pH)
CN <- CN %>% filter(Project.Year=="2021"&Horizon=="Organic") %>% dplyr::select(Sample.Date,Area,SiteName,PlotNumber,Rep,Horizon,N.Average,C.Average,CN.ratio) %>% mutate(Horizon = recode(Horizon, "Organic" = "O"))
names(CN) <- c("collection.date","Area","Site","Plot","Rep","Horizon","N%","C%","CN")

weights <- merge(wwdw, LOI, by=intersect(names(wwdw[,1:7]),names(LOI[,1:7])))

## remove analytical reps ###
weights_stacked <-gather(weights, "condition", "Value", 8:10)
weights_summary<-weights_stacked %>%
  group_by(collection.date, Area, Site,Plot,condition)%>%
  summarize(mean=mean(Value, na.rm = TRUE),
            n=length(Value),
            sd=sd(Value, na.rm = TRUE),
            se=(sd(Value, na.rm = TRUE)/sqrt(length(Value)) ),
            rsd=(sd(Value, na.rm = TRUE)/mean(Value))*100
  )

weights_means<-weights_summary %>%
  dplyr::select(collection.date, Area, Site,Plot,condition,mean) %>%
  spread(condition, mean)

CN_stacked <- gather(CN, "element","Value",7:9)
CN_summary<-CN_stacked %>%
  group_by(collection.date, Area, Site,Plot,element)%>%
  summarize(mean=mean(Value, na.rm = TRUE),
            n=length(Value),
            sd=sd(Value, na.rm = TRUE),
            se=(sd(Value, na.rm = TRUE)/sqrt(length(Value)) ),
            rsd=(sd(Value, na.rm = TRUE)/mean(Value))*100
  )

CN_means <- CN_summary %>% dplyr::select(collection.date,Area,Site,Plot,element,mean) %>%
  spread(element,mean)

soils_means_pt1 <- merge(weights_means, pH, by = intersect(names(weights_means[,1:4]), names(pH[,c(1:4)])))
soils_means <- merge(soils_means_pt1,CN_means, by = intersect(names(soils_means_pt1[,1:4]), names(CN_means[,1:4])))

soils_means<- soils_means %>% mutate(SampleEvent = recode(collection.date, "2-Jul-21"="Early Season", "5-Jul-21"="Early Season","19-Jul-21"="Mid Season","22-Jul-21"="Mid Season","2-Aug-21"="Late Season","5-Aug-21"="Late Season")) %>% dplyr::select(Area,Site,Plot,SampleEvent,pH,dw.ww,LOI,Moisture,`C%`,`N%`,CN) %>% mutate(Area=recode(Area, "East"="Acidic","West"="Non-acidic"))
soils_means[,6:11] <- round(soils_means[,6:11], 2)
names(soils_means) <- c("Area","SiteName","PlotNumber","SampleEvent","pH","dw.ww","LOI","Moisture","C%","N%","CN")
write.csv(soils_means, "formatted spreadsheets/soil_conditions_plot.csv",row.names = F)


# thaw depth data #
thaw <- read.csv("raw data/PhIr2021_thaw.csv")
thaw_stacked <- gather(thaw, "Obs", "Value", 6:20)

# filter data to dates when we collected soil cores
thaw_sample <- thaw_stacked %>% filter(Date == "1-Jul-2021" | Date == "24-Jul-2021" | Date == "7-Aug-2021") %>% mutate(SampleEvent = recode(Date, "1-Jul-2021"="Early Season","24-Jul-2021"="Mid Season","7-Aug-2021"="Late Season")) %>% dplyr::select(SampleEvent,Area,Site,Value) %>% mutate(Area=recode(Area, "East"="Acidic","West"="Non-acidic"))

write.csv(thaw_sample, "formatted spreadsheets/thaw_sample.csv", row.names = F)

### summarize PSI data ####
PSI <- read.csv("raw data/PhIr2021_PSI_2.csv", na.strings=c("#DIV/0!", "#VALUE!", " ", "NA"))
names(PSI) <- c("Project", "SoilSampleDate", "Area", "SiteName", "PlotNumber", "SampleRep", "Horizon", "SampleName", "PSI")
PSI <- PSI %>% filter(!is.na(PlotNumber))%>% mutate(Area=recode(Area, "West"="Non-acidic","East"="Acidic")) 

PSI_Summary<-PSI %>%
  group_by(SoilSampleDate, Area, SiteName, PlotNumber, Horizon)%>%
  summarize(mean=mean(PSI, na.rm = TRUE),
            n=length(PSI),
            sd=sd(PSI, na.rm = TRUE),
            se=(sd(PSI, na.rm = TRUE)/sqrt(length(PSI)) ),
            rsd=(sd(PSI, na.rm = TRUE)/mean(PSI))*100
  )

 PSI_Summary <- PSI_Summary %>% mutate(SampleEvent = recode(SoilSampleDate, "2-Jul-21"="Early Season", "5-Jul-21"="Early Season","19-Jul-21"="Mid Season","22-Jul-21"="Mid Season","2-Aug-21"="Late Season","5-Aug-21"="Late Season"))

 
### compile spreadsheets into master spreadsheet #### 
compiled_pfrac <- Pfrac_Means %>% dplyr::select(SampleEvent,Area,SiteName,PlotNumber,Horizon,PfracBDFe,PfracBDNRP,PfracBDSRP,PfracH2ONRP,PfracH2OSRP,PfracHClCa,PfracHClP,PfracNaOHAl,PfracNaOHNRP,PfracNaOHSRP,PfracNRP,PfracRP,PfracTotAl,PfracTotCa,PfracTotFe,PfracTotP,ResP)
 compiled_PSI <- PSI_Summary %>% dplyr::select(SampleEvent,Area,SiteName,PlotNumber,Horizon,mean)
 compiled_fefrac<- Fefrac_Means %>% dplyr::select(SampleEvent,Area,SiteName,PlotNumber,Horizon,FefracDHFe,FefracHHFe,FefracPPFe,ResFe)
 compiled_soils <- soils_means %>% dplyr::select(SampleEvent,Area,SiteName,PlotNumber,LOI,Moisture,pH,`C%`,`N%`,CN)
 
 compiled_pt1 <- merge(compiled_pfrac,compiled_PSI, by = intersect(names(compiled_pfrac[,1:6]),names(compiled_PSI[,1:6])))
compiled_pt2 <- merge(compiled_fefrac,compiled_soils,by=intersect(names(compiled_fefrac[,2:5]),names(compiled_soils[,1:4]))) 
compiled_pt3 <- merge(compiled_pt1,compiled_pt2, by = intersect(names(compiled_pt1[,1:6]),names(compiled_pt2[,c(5,1,2,3,4,6)])))

soils_PSI <- compiled_pt3 %>% ungroup() %>% mutate(treatment=paste(Area,SiteName, sep = " ")) %>% dplyr::select(SampleEvent,Area,SiteName,treatment,PlotNumber,Horizon,mean,pH,LOI,Moisture,`C%`,`N%`,CN,PfracBDFe, PfracBDNRP,PfracBDSRP,PfracH2ONRP,PfracH2OSRP, PfracHClCa,PfracHClP, PfracNaOHAl,PfracNaOHNRP,PfracNaOHSRP,PfracNRP,PfracRP,PfracTotAl,PfracTotCa,PfracTotFe, PfracTotP,ResP,FefracDHFe,FefracHHFe,FefracPPFe,ResFe)
names(soils_PSI) <- c("SampleEvent","Area", "SiteName", "treatment","PlotNumber","Horizon","PfracPSI","pH","LOI","moisture","C%","N%","CN_ratio","PfracBDFe", "PfracBDNRP","PfracBDSRP","PfracH2ONRP","PfracH2OSRP", "PfracHClCa","PfracHClP", "PfracNaOHAl","PfracNaOHNRP","PfracNaOHSRP","PfracNRP","PfracRP","PfracTotAl","PfracTotCa","PfracTotFe", "PfracTotP","ResP","FefracDHFe","FefracHHFe","FefracPPFe","ResFe")

write.csv(soils_PSI,"formatted spreadsheets/soils_PSI.csv", row.names = F)
