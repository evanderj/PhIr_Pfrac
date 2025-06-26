setwd("~/PhIr/PhIr 2021/For github")
source("0-packages.R")

soils_PSI <- read.csv("formatted spreadsheets/soils_PSI.csv")
rP_stacked <- read.csv("formatted spreadsheets/Pfrac_RP_stacked.csv")
TP_stacked <- read.csv("formatted spreadsheets/Pfrac_TP_stacked.csv")
iron_stacked <- read.csv("formatted spreadsheets/Fefrac_Fe_stacked.csv")
thaw <- read.csv("formatted spreadsheets/thaw_sample.csv")


soils_PSI <- soils_PSI %>% mutate(SiteName= factor(SiteName, levels = c("Dry", "Mesic", "Hydric"))) %>% mutate(Area=factor(Area)) %>% mutate(SampleEvent=factor(SampleEvent,levels = c("Early Season", "Mid Season", "Late Season"))) %>% mutate(treatment=factor(treatment,levels=c("Acidic Dry", "Acidic Mesic", "Acidic Hydric", "Non-acidic Dry", "Non-acidic Mesic", "Non-acidic Hydric")))
soils_PSI$lnH2OSRP <- log(soils_PSI$PfracH2OSRP)

rP_stacked <- rP_stacked %>% mutate(SiteName= factor(SiteName, levels = c("Dry", "Mesic", "Hydric"))) %>% mutate(Area=factor(Area)) %>% mutate(Fraction=factor(Fraction, levels = c("PfracH2OSRP","PfracBDSRP","PfracNaOHSRP","PfracHClP")))
TP_stacked <- TP_stacked %>% mutate(SiteName= factor(SiteName, levels = c("Dry", "Mesic", "Hydric"))) %>% mutate(Area=factor(Area)) %>% mutate(Fraction=factor(Fraction, levels = c("PfracRP","PfracNRP","ResP")))
iron_stacked <- iron_stacked %>% mutate(SiteName= factor(SiteName, levels = c("Dry", "Mesic", "Hydric"))) %>% mutate(Area=factor(Area)) %>% mutate(Fraction=factor(Fraction, levels = c("FefracPPFe","FefracHHFe","FefracDHFe")))

thaw<- thaw %>% mutate(Site= factor(Site, levels = c("Dry", "Mesic", "Hydric")))%>% mutate(Area=factor(Area)) %>% mutate(SampleEvent=factor(SampleEvent, levels = c("Early Season", "Mid Season", "Late Season")))


## soil conditions supplemental ##
thaw_box <- thaw %>% mutate(Site = recode(Site, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))%>% mutate(Site=factor(Site, levels = c("Upland","Midslope","Lowland"))) %>% ggplot(aes(x = Site, y = Value)) +
  geom_boxplot(aes(fill = SampleEvent), show.legend = F)+
  scale_fill_manual(values = pnw_palette(name = "Sunset", 3))+
  facet_grid(cols = vars(Area))+
  scale_y_reverse()+
  theme_er2()+
  xlab("") + ylab("Thaw Depth (cm)")

moist_summ <- soils_PSI %>%  group_by(SampleEvent, Area, SiteName, treatment) %>% 
  summarize(mean=mean(moisture, na.rm = TRUE),
            sd=sd(moisture, na.rm = TRUE),
            se=(sd(moisture, na.rm = TRUE)/sqrt(length(moisture)) ),
            rsd=(sd(moisture, na.rm = TRUE)/mean(moisture))*100)

moisture <- moist_summ %>% mutate(SiteName = recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland")) %>% ggplot(aes(x = SiteName, ymin = mean-sd, ymax=mean+sd,fill = SampleEvent))+
  geom_bar(position = "dodge", stat = "identity", aes(y = mean))+
  geom_errorbar(position = position_dodge(width = 0.9), width = 0.4)+
  scale_fill_manual(values = PNWColors::pnw_palette("Sunset", 3))+
  coord_cartesian(ylim= c(55,90))+
  facet_grid(.~Area)+
  ylab("Moisture (%)")+
  xlab("Hillslope Position")+
  theme_er2()+
  theme(legend.position = "bottom")

legend <- get_legend(moisture)
moisture <- moisture + theme(legend.position="none")
grid.arrange(legend, thaw_box, moisture, heights = c(2,5,5))


### total soil phosphorus ##
TP_stacked %>% mutate(SiteName = recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland")) %>%
  ggplot(aes(x=SiteName,y=mean, fill = Fraction)) +
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(aes(ymax = Accumulation + se, ymin = Accumulation - se),width = 0.5, position = position_dodge2(width = 0.3))+
  theme_er2()+
  facet_grid(cols = vars(Area)) + 
  xlab("Hillslope Position") + ylab("Phosphorus (µg P/g-soil)") + 
  scale_fill_manual(values = PNWColors::pnw_palette("Shuksan",5), labels = c("Reactive P", "Non-reactive P", "Residual P"))
ggsave("figures/TP_stacked_snap.png", width=10, height=5, units=c("in"))

### reactive P ##
rP_stacked %>% mutate(SiteName=recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))%>%
  ggplot(aes(x=SiteName,y=mean, fill = Fraction)) + 
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(aes(ymax = Accumulation + se, ymin = Accumulation - se), width = 0.5, position = position_dodge(width = 0.3)) +
  facet_grid(cols = vars(Area)) +
  theme_er2()+
  scale_fill_manual(values = c("#fbdfa2","#de9b71","#41476b","#9e6374"), labels = c("Loosely sorbed rP", "rP~iron oxides", "rP~aluminum oxides", "P~calcareous minerals"))  + 
  labs(x = "Hillslope Position", y = "reactive P (µg P/g-soil)")
ggsave("figures/RP_stacked_snap.png", width=10, height=5, units=c("in"))


### iron ###
iron_stacked %>% mutate(SiteName = recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))%>%
  ggplot(aes(x=SiteName,y=mean, fill = Fraction)) + 
  geom_bar(position = "stack", stat = "identity") +
  geom_errorbar(aes(ymax = Accumulation + se, ymin = Accumulation - se), width = 0.5, position = position_dodge2(width = 0.3)) +
  scale_fill_manual(values = c("#ecc0a1","#c26a7a","#984136"), labels = c("iron~soil organic matter","poorly crystalline iron","crystalline iron")) +
  theme_er2()+
  labs(x = "Hillslope Position", y = "Fe (mg Fe/g-soil)") + 
  facet_grid(cols = vars(Area))
ggsave("figures/Fe_stacked_snap.png", width=10, height=5, units=c("in"))



### P reationships ###
SRP_H2O_D <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI[which(soils_PSI$SiteName=="Dry"),])
base::summary(SRP_H2O_D)
start_D <- list(a = exp(coef(SRP_H2O_D)[1]), b = coef(SRP_H2O_D)[2])

SRP_H2O_M <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI[which(soils_PSI$SiteName=="Mesic"),])
base::summary(SRP_H2O_M)
start_M <- list(a = exp(coef(SRP_H2O_M)[1]), b = coef(SRP_H2O_M)[2])

SRP_H2O_H <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI[which(soils_PSI$SiteName=="Hydric"),])
base::summary(SRP_H2O_H)
start_H <- list(a = exp(coef(SRP_H2O_H)[1]), b = coef(SRP_H2O_H)[2])

nls(PfracH2OSRP~(a*exp(b*PfracBDSRP)), data = soils_PSI[which(soils_PSI$SiteName == "Dry"),], start = start_D)
nls(PfracH2OSRP~(a*exp(b*PfracBDSRP)), data = soils_PSI[which(soils_PSI$SiteName == "Mesic"),], start = start_M)
nls(PfracH2OSRP~(a*exp(b*PfracBDSRP)), data = soils_PSI[which(soils_PSI$SiteName == "Hydric"),], start = start_H)

best_fit <- nls(PfracH2OSRP~(a*exp(b*PfracBDSRP)), data = soils_PSI, start = start_D)
best_fit2 <- lm(lnH2OSRP~PfracBDSRP, data = soils_PSI)

exp_fun <- data.frame(SiteName = c("Upland","Midslope","Lowland","Upland","Midslope","Lowland"), X = c(20,20,20,15,15,15), Y = c(30,30,30,25,25,25), lab = c( "y == 0.195*e^{0.133*x}","y == 1.99*e^{0.046*x}","y == 0.416*e^{0.047*x}","p == 0.003","p == 0.004","p < 0.001"))
exp_fun$SiteName <- factor(exp_fun$SiteName, levels = c("Upland","Midslope","Lowland"))

soils_PSI %>% mutate(SiteName = recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))%>%
  ggplot(aes(x=PfracBDSRP, y=PfracH2OSRP))+
  geom_point(aes(color = Area), size = 3)+
  geom_text(data = exp_fun, aes(x = X, y = Y,label = lab), parse = TRUE)+
  geom_smooth(method = "nls", formula = y~a*exp(b*x), se = FALSE, method.args = list(start=coef(best_fit)), color="#59629b", linewidth =0.8) +
  xlab("rP~iron oxides (µg P/g-soil)")+
  ylab("loosely sorbed rP (µg P/g-soil)")+
  scale_color_manual(values=pnw_palette("Starfish",2), name = "soil pH") + 
  theme_er2()+
  facet_grid(.~SiteName)
ggsave("figures/H2O~BDSRP_fit.png", width = 9, height = 5, units = "in")


### H2OSRP and iron ##
H2O_iron <- soils_PSI %>% dplyr::select(Area,SiteName,PfracH2OSRP,PfracBDFe,FefracDHFe,FefracHHFe,FefracPPFe) %>% gather("Fraction", "Value", c(4:7)) %>% mutate(Fraction=factor(Fraction, levels = c("PfracBDFe", "FefracPPFe","FefracHHFe","FefracDHFe"), label = c("iron oxides", "iron~soil organic matter","poorly \ncrystalline iron","crystalline iron")))

H2O_iron %>% mutate(SiteName=recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))%>%
  ggplot(aes(x = Value, y = PfracH2OSRP, color = SiteName,shape = Area))+
  geom_point(size = 4)+
  scale_color_manual(values = c("#7F6552","#6A9741","#5A8093"), name = "Hillslope\nPosition")+
  scale_shape(name = "soil pH")+
  facet_grid(.~Fraction, scales = "free")+
  labs(y= "loosely sorbed rP (µg P/g-soil)", x = "Iron concentration (mg Fe/g-soil)")+
  facet_wrap(.~Fraction, scales = "free")+
  theme_er2()+
  theme(strip.text.x = element_text(size=16))
ggsave("figures/iron~H2O.png", width = 12, height = 7, units = "in")


### BDSRP and iron ##
Fe_HH <- soils_PSI %>% mutate(SiteName = recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland")) %>% ggplot(aes(x=FefracHHFe, y=PfracBDSRP)) +
  geom_point(size = 4, aes(color = SiteName, shape = Area), show.legend = F)+
  scale_color_manual(values = c("#7F6552","#6A9741","#5A8093"), name = "Hillslope\nPosition")+
  scale_shape_manual(values = c(16,17),name = "soil pH")+
  labs(x = "poorly crystalline iron\n (mg Fe/g-soil)", y = "") +
  theme_er2()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))

Fe_PP<- soils_PSI %>% mutate(SiteName = recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland")) %>% ggplot(aes(x=FefracPPFe, y=PfracBDSRP)) +
  geom_point(size = 4, aes(color= SiteName, shape = Area))+
  scale_color_manual(values = c("#7F6552","#6A9741","#5A8093"), name = "Hillslope\nPosition")+
  scale_shape_manual(values = c(16,17),name = "soil pH")+
  labs(x = "iron~soil organic matter\n (mg Fe/g-soil)", y ="") +
  theme_er2()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))


legend <- get_legend(Fe_PP)
Fe_PP <- Fe_PP + theme(legend.position="none")

grid.arrange(Fe_PP, Fe_HH, legend, ncol = 3, widths = c(6,6,2), left = textGrob("rP~iron oxides (µg P/g-soil)", gp=gpar(fontsize=20, font=3),rot = 90))


### SRP with aluminum and calcium ##
Al_SRP <-soils_PSI %>% filter(Area=="Acidic") %>% mutate(SiteName=recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))%>% ggplot(aes(x = PfracNaOHAl, y = PfracNaOHSRP, color = SiteName))+
  geom_point(size = 5, show.legend = F)+
  scale_color_manual(values = c("#7F6552","#6A9741","#5A8093"), name = "Hillslope position")+
  labs(y ="rP~Aluminum oxides (µg P/g-soil)", x ="Aluminum oxides (mg Al/g-soil)", subtitle = "Acidic tundra")+
  theme_er2()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))


Ca_P <- soils_PSI %>% filter(Area=="Non-acidic") %>% mutate(SiteName=recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))%>%  ggplot(aes(x = PfracHClCa, y = PfracHClP, color = SiteName))+
  geom_point(size = 5)+
  scale_color_manual(values = c("#7F6552","#6A9741","#5A8093"), name = "Hillslope \nPosition")+
  labs(y = "P~Calcareous minerals (µg P/g-soil)", x = "Calcium carbonate (mg Ca/g-soil)", subtitle = "Non-acidic tundra")+
  theme_er2()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))

legend <- get_legend(Ca_P)
Ca_P <- Ca_P + theme(legend.position = "none")

grid.arrange(Al_SRP, Ca_P,legend, widths = c(6,6,2))


### PSI and H2OSRP ##
soils_PSI %>% mutate(SiteName=recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))%>%
  ggplot(aes(x = PfracPSI, y = PfracH2OSRP, color = SiteName))+
  geom_point(size = 4)+
  scale_color_manual(values =c("#7F6552","#6A9741","#5A8093"), name = "Hillslope position")+
  labs( x= "Phosphate sorption index", y = "loosely sorbed rP (µg P/g-soil)")+
  theme_er2()
ggsave("figures/PSI~H2OrP.png", width = 10, height = 5, units = "in")


### PSI, iron and calcium ##
PSI_Ca <- soils_PSI %>% filter(Area == "Non-acidic") %>% mutate(SiteName=recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))%>% ggplot(aes(x = PfracHClCa, y = PfracPSI, color = SiteName))+
  geom_point(size = 5, show.legend = F)+
  scale_color_manual(values = c("#7F6552","#6A9741","#5A8093"), name = "Hillslope position")+
  labs(y = "", x = "Calcium carbonate (mg Ca/g-soil)", subtitle = "Non-acidic tundra")+
  scale_y_continuous(limits = c(0,140))+
  theme_er2()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))


PSI_Fe <-soils_PSI %>% filter(Area == "Acidic") %>% mutate(SiteName=recode(SiteName, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))%>% ggplot(aes(x = PfracBDFe, y = PfracPSI, color = SiteName))+
  geom_point(size = 5)+
  scale_color_manual(values = c("#7F6552","#6A9741","#5A8093"), name = "Hillslope \nposition")+
  scale_y_continuous(limits = c(0,140))+
  labs(x = "Iron oxides (mg Fe/g-soil)", y = "",subtitle = "Acidic tundra")+
  theme_er2()+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))


legend <- get_legend(PSI_Fe)
PSI_Fe <- PSI_Fe + theme(legend.position = "none")

grid.arrange(PSI_Fe,PSI_Ca,legend, ncol=3,widths=c(6,6,2),left = textGrob("Phosphate sorption index", gp=gpar(fontsize=20, font=3),rot = 90))
