setwd("~/PhIr/PhIr 2021")
source("scripts/0-packages.R")

moisture<-read.csv("PhIr 2021 soil conditions/PhIr_redox_porechem-main/processed/PhIr_moisture_temp_salinity_2021_2022.csv")
redox<-read.csv("PhIr 2021 soil conditions/PhIr_redox_porechem-main/processed/PhIr_redox_2021_2022.csv")

## fix timestamps and sample IDs ##
moisture$datetime <- strptime(moisture$TIMESTAMP, format = "%m/%d/%Y %H:%M")
moisture_fixed <- moisture %>% filter(datetime < "2022-01-01 00:00:00")%>% filter(moisture > 10&temp>0) %>%mutate(area= recode(area, "east"="Acidic", "west"= "Non-acidic")) %>% mutate(site=recode(site, "hydric"="Hydric","mesic"="Mesic","dry"="Dry")) %>% mutate(site=factor(site, levels = c("Dry", "Mesic", "Hydric"))) %>% mutate(treatment=paste(area,site)) %>% mutate(sample=if_else(area=="Acidic"&datetime>="2021-06-25 00:00:00"&datetime<="2021-07-02 07:00:00","early",ifelse(area=="Acidic"&datetime>="2021-07-12 00:00:00"&datetime<="2021-07-19 07:00:00", "mid",ifelse(area=="Acidic"&datetime>="2021-07-26 00:00:00"&datetime<="2021-08-02 07:00:00","late",ifelse(area=="Non-acidic"&datetime>="2021-06-28 00:00:00"&datetime<="2021-07-05 07:00:00","early",ifelse(area=="Non-acidic"&datetime>="2021-07-15 00:00:00"&datetime<="2021-07-22 07:00:00","mid",ifelse(area=="Non-acidic"&datetime>="2021-07-29 00:00:00"&datetime<="2021-08-05 07:00:00","late",NA))))))) %>% mutate(sample=factor(sample,levels = c("early","mid","late",NA)))


redox_fixed <- redox %>% mutate(datetime= strptime(redox$TIMESTAMP, format = "%m/%d/%Y %H:%M")) %>% filter(datetime < "2022-01-01"&condition =="unfrozen"&redox_mV >-250&redox_mV<750)%>% mutate(area=recode(area, "acidic tundra"="Acidic","non-acidic tundra"="Non-acidic")) %>% mutate(site=recode(site, "hydric"="Hydric","mesic"="Mesic","dry"="Dry")) %>% mutate(site=factor(site,levels = c("Dry","Mesic","Hydric"))) %>% mutate(treatment=paste(area,site)) %>% mutate(treatment=factor(treatment, levels = c("Acidic Dry", "Acidic Mesic", "Acidic Hydric", "Non-acidic Dry", "Non-acidic Mesic", "Non-acidic Hydric"))) %>% mutate(sample=if_else(area=="Acidic"&datetime>="2021-06-25 00:00:00"&datetime<="2021-07-02 07:00:00","early",ifelse(area=="Acidic"&datetime>="2021-07-12 00:00:00"&datetime<="2021-07-19 07:00:00", "mid",ifelse(area=="Acidic"&datetime>="2021-07-26 00:00:00"&datetime<="2021-08-02 07:00:00","late",ifelse(area=="Non-acidic"&datetime>="2021-06-28 00:00:00"&datetime<="2021-07-05 07:00:00","early",ifelse(area=="Non-acidic"&datetime>="2021-07-15 00:00:00"&datetime<="2021-07-22 07:00:00","mid",ifelse(area=="Non-acidic"&datetime>="2021-07-29 00:00:00"&datetime<="2021-08-05 07:00:00","late",NA))))))) %>% mutate(sample=factor(sample,levels = c("early","mid","late",NA)))


moisture_clipped<- moisture_fixed[which(moisture_fixed$datetime >= "2021-07-01 00:00:00" & moisture_fixed$datetime<= "2021-08-06 00:00:00"),]
redox_clipped<- redox_fixed[which(redox_fixed$datetime >= "2021-07-01 00:00:00" & redox_fixed$datetime <= "2021-08-06 00:00:00"),]


## moisture_redox composite figure ####
#averaging moisture across depths
moisture_summary<-moisture_clipped %>%  filter(depth_cm < 20)%>%
  group_by(datetime, treatment)%>%
  summarize(mean=mean(moisture, na.rm = TRUE),
            n=length(moisture),
            sd=sd(moisture, na.rm = TRUE),
            se=(sd(moisture, na.rm = TRUE)/sqrt(n())),
            rsd=(sd(moisture, na.rm = TRUE)/mean(moisture, na.rm = T))*100
  )

# filling gaps in continuous data with NA's
moisture_av_time <- moisture_summary[which(moisture_summary$treatment != "Non-acidic Hydric"),1:3] %>% pivot_wider(names_from = treatment, values_from = mean)
moisture_av_time$datetime <- as.character(moisture_av_time$datetime)
NA_fill <- data.frame(datetime = seq.POSIXt(from= as.POSIXct(moisture_av_time$datetime[1]), to = as.POSIXct(moisture_av_time$datetime[length(moisture_av_time$datetime)]), by = "15 mins"))
NA_fill$datetime <- as.character(NA_fill$datetime)
NA_merge <- merge(moisture_av_time, NA_fill, by.y = "datetime", all.y= TRUE)
moisture_se_time <- moisture_summary[which(moisture_summary$treatment != "Non-acidic Hydric"),c(1,2,6)] %>% pivot_wider(names_from = treatment, values_from = se)
moisture_se_time$datetime <- as.character(moisture_se_time$datetime)
NA_se_merge<-merge(moisture_se_time,NA_fill, by.y = "datetime",all.y =TRUE)
moisture_av_fill <- NA_merge %>% gather("treatment", "mean", 2:6)
moisture_se_fill <- NA_se_merge %>% gather("treatment","se",2:6)
moisture_fill <- merge(moisture_av_fill,moisture_se_fill, by = intersect(names(moisture_av_fill[,1:2]),names(moisture_se_fill[,1:2])))

nah_av_time <- moisture_summary[which(moisture_summary$treatment == "Non-acidic Hydric"),1:3] %>% pivot_wider(names_from = treatment, values_from = mean)
nah_av_time$datetime<-as.character(nah_av_time$datetime)
nah_fill <- data.frame(datetime = seq.POSIXt(from= as.POSIXct(nah_av_time$datetime[1]), to = as.POSIXct(nah_av_time$datetime[length(nah_av_time$datetime)]), by = "15 mins"))
nah_fill$datetime <- as.character(nah_fill$datetime)
nah_merge <- merge(nah_av_time, nah_fill, by.y = "datetime", all.y = TRUE)
moisture_av_nah <- nah_merge %>% gather("treatment", "mean", 2)

nah_se_time <-   moisture_summary[which(moisture_summary$treatment == "Non-acidic Hydric"),c(1,2,6)] %>% pivot_wider(names_from = treatment, values_from = se)
nah_se_time$datetime <- as.character(nah_se_time$datetime)
nah_se_merge <- merge(nah_se_time,nah_fill, by.y="datetime",all.y=TRUE)
moisture_se_nah <- nah_se_merge %>% gather("treatment","se",2)
moisture_nah <- merge(moisture_av_nah, moisture_se_nah, by=intersect(names(moisture_av_nah[,1:2]),names(moisture_se_nah[,1:2])))

moisture_filled <- rbind(moisture_fill, moisture_nah)
moisture_filled <- moisture_filled %>% mutate(area = if_else(treatment %in% c("Acidic Dry", "Acidic Mesic","Acidic Hydric"), "Acidic", "Non-acidic"))  %>% mutate(site=if_else(treatment %in% c("Acidic Dry", "Non-acidic Dry"), "Dry",ifelse(treatment %in% c("Acidic Mesic","Non-acidic Mesic"), "Mesic","Hydric"))) %>% mutate(treatment=factor(treatment,levels = c("Acidic Dry", "Acidic Mesic", "Acidic Hydric", "Non-acidic Dry","Non-acidic Mesic", "Non-acidic Hydric"))) %>% mutate(site=factor(site, levels=c("Dry","Mesic","Hydric")))

# moisture figure
moisture_comp <-ggplot(moisture_filled)+
  geom_line(aes(y = mean, x = as.POSIXct(strptime(datetime, format="%Y-%m-%d %H:%M:%S")), color = site), linewidth = 1, show.legend = FALSE)+
  geom_ribbon(aes(x = as.POSIXct(strptime(datetime, format="%Y-%m-%d %H:%M:%S")),y = mean, ymin = mean-se, ymax = mean+se, fill = site), alpha = 0.3,show.legend = FALSE)+
  scale_x_datetime(date_labels = "%b-%d")+
  scale_color_manual(values =c("#7F6552","#6A9741","#5A8093"), name = "soil condition")+
  labs(x = " ", y = "moisture %")+
  scale_fill_manual(values =c("#7F6552","#6A9741","#5A8093"), name = "soil condition")+
  geom_vline(data = subset(moisture_filled, area == "Acidic"), aes(xintercept = as.POSIXct("2021-07-02 07:30:00")),linetype="dashed",color="black")+
  geom_vline(data = subset(moisture_filled, area == "Acidic"), aes(xintercept = as.POSIXct("2021-07-19 07:30:00")),linetype="dashed",color="black")+
  geom_vline(data = subset(moisture_filled, area == "Acidic"), aes(xintercept = as.POSIXct("2021-08-02 07:30:00")),linetype="dashed",color="black")+
  geom_vline(data = subset(moisture_filled, area == "Non-acidic"), aes(xintercept = as.POSIXct("2021-07-05 07:30:00")),linetype="dashed",color="black")+
  geom_vline(data = subset(moisture_filled, area == "Non-acidic"), aes(xintercept = as.POSIXct("2021-07-21 07:30:00")),linetype="dashed",color="black")+
  geom_vline(data = subset(moisture_filled, area == "Non-acidic"), aes(xintercept = as.POSIXct("2021-08-05 07:30:00")),linetype="dashed",color="black")+
  facet_grid(.~area)+
  theme_er2()


# averaging redox across depths
redox_summary<-redox_clipped %>% filter(depth_cm >= 5 & depth_cm <= 15)%>%
  group_by(datetime, treatment)%>%
  summarize(mean=mean(redox_mV, na.rm = TRUE),
            n=length(redox_mV),
            sd=sd(redox_mV, na.rm = TRUE),
            se=(sd(redox_mV, na.rm = TRUE)/sqrt(n())),
            rsd=(sd(redox_mV, na.rm = TRUE)/mean(redox_mV))*100)

# filling gaps in continuous data with NA's
redox_av_time <- redox_summary[which(redox_summary$treatment != "Non-acidic Hydric"),1:3] %>% pivot_wider(names_from = treatment, values_from = mean)
NA_fill <- data.frame(datetime = seq.POSIXt(from= as.POSIXct(redox_av_time$datetime[1]), to = as.POSIXct(redox_av_time$datetime[length(redox_av_time$datetime)]), by = "15 mins"))
redox_av_time$datetime <- as.POSIXct(strptime(redox_av_time$datetime, format = "%Y-%m-%d %H:%M:%S"))
NA_merge <- merge(redox_av_time, NA_fill, by.y = "datetime", all.y= TRUE)
redox_av_fill <- NA_merge %>% gather("treatment", "mean", 2:6)

redox_se_time <- redox_summary[which(redox_summary$treatment != "Non-acidic Hydric"),c(1,2,6)] %>% pivot_wider(names_from = treatment, values_from = se)
redox_se_time$datetime <-as.POSIXct(strptime(redox_se_time$datetime, format = "%Y-%m-%d %H:%M:%S"))
NA_se_merge <- merge(redox_se_time, NA_fill, by.y = "datetime", all.y= TRUE)
redox_se_fill <- NA_se_merge %>% gather("treatment", "se", 2:6)
redox_fill <-merge(redox_av_fill,redox_se_fill,by=intersect(names(redox_av_fill[,1:2]),names(redox_se_fill[,1:2])))


nah_av_time <-   redox_summary[which(redox_summary$treatment == "Non-acidic Hydric"),1:3] %>% pivot_wider(names_from = treatment, values_from = mean)
nah_fill <- data.frame(datetime = seq.POSIXt(from= as.POSIXct(nah_av_time$datetime[1]), to = as.POSIXct(nah_av_time$datetime[length(nah_av_time$datetime)]), by = "15 mins"))
nah_av_time$datetime <- as.POSIXct(strptime(nah_av_time$datetime, format = "%Y-%m-%d %H:%M:%S"))
nah_merge <- merge(nah_av_time, nah_fill, by.y = "datetime", all.y = TRUE)
redox_av_nah <- nah_merge %>% gather("treatment", "mean", 2)

nah_se_time <-   redox_summary[which(redox_summary$treatment == "Non-acidic Hydric"),c(1,2,6)] %>% pivot_wider(names_from = treatment, values_from = se)
nah_se_time$datetime <- as.POSIXct(strptime(nah_se_time$datetime, format = "%Y-%m-%d %H:%M:%S"))
nah_se_merge <- merge(nah_se_time, nah_fill, by.y = "datetime", all.y = TRUE)
redox_se_nah <- nah_se_merge %>% gather("treatment", "se", 2)
redox_nah <- merge(redox_av_nah,redox_se_nah, by=intersect(names(redox_av_nah[,1:2]),names(redox_se_nah[,1:2])))

redox_filled <- rbind(redox_fill, redox_nah)
redox_filled <- redox_filled %>% mutate(area = if_else(treatment %in% c("Acidic Dry", "Acidic Mesic","Acidic Hydric"), "Acidic", "Non-acidic")) %>% mutate(site=if_else(treatment %in% c("Acidic Dry", "Non-acidic Dry"), "Dry",ifelse(treatment %in% c("Acidic Mesic","Non-acidic Mesic"), "Mesic","Hydric"))) %>% mutate(treatment=factor(treatment,levels = c("Acidic Dry", "Acidic Mesic", "Acidic Hydric", "Non-acidic Dry","Non-acidic Mesic", "Non-acidic Hydric"))) %>% mutate(site=factor(site, levels=c("Dry","Mesic","Hydric"))) %>% mutate(site = recode(site, "Dry"="Upland","Mesic"="Midslope","Hydric"="Lowland"))

# redox figure
redox_comp <-redox_filled %>% ggplot()+
  geom_line(aes(y = mean, x = as.POSIXct(datetime), color = site), linewidth = 1)+
  geom_ribbon(aes(x = as.POSIXct(datetime),y = mean, ymin = mean-se, ymax = mean+se, fill = site), alpha = 0.3)+
  scale_x_datetime(date_labels = "%b-%d")+
  scale_color_manual(values =c("#7F6552","#6A9741","#5A8093"), name ="hillslope \nposition")+
  scale_fill_manual(values =c("#7F6552","#6A9741","#5A8093"), name = "hillslope \nposition")+
  labs(x = " ", y = "Redox potential (mV)")+
  geom_vline(data= subset(redox_filled, area == "Acidic"), aes(xintercept = as.POSIXct("2021-07-02 07:30:00")),linetype="dashed",color="black")+
  geom_vline(data= subset(redox_filled, area == "Acidic"), aes(xintercept = as.POSIXct("2021-07-19 07:30:00")),linetype="dashed",color="black")+
  geom_vline(data= subset(redox_filled, area == "Acidic"), aes(xintercept = as.POSIXct("2021-08-02 07:30:00")),linetype="dashed",color="black")+
  geom_vline(data= subset(redox_filled, area == "Non-acidic"), aes(xintercept = as.POSIXct("2021-07-05 07:30:00")),linetype="dashed",color="black")+
  geom_vline(data= subset(redox_filled, area == "Non-acidic"), aes(xintercept = as.POSIXct("2021-07-21 07:30:00")),linetype="dashed",color="black")+
  geom_vline(data= subset(redox_filled, area == "Non-acidic"), aes(xintercept = as.POSIXct("2021-08-05 07:30:00")),linetype="dashed",color="black")+
  facet_grid(.~area)+
  theme_er2()

legend <- get_legend(redox_comp)
redox_comp <- redox_comp + theme(legend.position="none")

grid.arrange(arrangeGrob(moisture_comp, redox_comp), legend, ncol = 2, widths = c(7,1))


### Manuscript values ####
redox_summ <- redox_clipped %>% filter(depth_cm >=5 &depth_cm<=15) %>% 
  group_by(treatment) %>%
  reframe(mean= mean(redox_mV, na.rm=T),
          se= sd(redox_mV, na.rm = T)/sqrt(n()),
          n = length(redox_mV))

print(redox_summ)

redox_A_treat <- lm(redox_mV~treatment, data = redox_clipped[redox_clipped$area=="Acidic"&redox_clipped$depth_cm<=15&redox_clipped$depth_cm>=5,])
summary(redox_A_treat)

redox_NA_treat <- lm(redox_mV~treatment, data = redox_clipped[redox_clipped$area=="Non-acidic"&redox_clipped$depth_cm<=15&redox_clipped$depth_cm>=5,])
summary(redox_NA_treat)


summary(redox_filled$mean[which(redox_filled$area =="Acidic"& redox_filled$datetime >= as.POSIXct("2021-07-01 00:00:00") & redox_filled$datetime <= as.POSIXct("2021-07-08 00:00:00"))],na.rm = T)
summary(redox_filled$mean[which(redox_filled$area =="Acidic"& redox_filled$datetime >= as.POSIXct("2021-07-15 00:00:00") & redox_filled$datetime <= as.POSIXct("2021-07-20 00:00:00"))],na.rm = T)


