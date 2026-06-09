#packages for data vis and manipulating data
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggthemes)
library(stats)
library(base)
library(PNWColors)
library(NatParksPalettes)
library(grid)
library(gridExtra)
library(metR)
library(HH)
library(AICcmodavg)
library(pls)
library(vegan)
library(ggfortify)
library(leaps)
library(broom)


theme_er2 <- function() {  # this for all the elements common across plots
  theme_bw() %+replace%
    theme(legend.key=element_blank(),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15),
          legend.key.size = unit(1.5, 'lines'),
          panel.border = element_rect(color="black",linewidth =2, fill = NA),
          
          plot.title = element_text(hjust = 0.5, size = 21),
          plot.subtitle = element_text(hjust = 0.5, size = 17, lineheight = 1.5),
          axis.text = element_text(size = 14, color = "black"),
          axis.title = element_text(size = 14, face = "bold", color = "black"),
          
          # formatting for facets
          panel.background = element_blank(),
          strip.background = element_rect(colour="white", fill="white"), #facet formatting
          panel.spacing.x = unit(1.5, "lines"), #facet spacing for x axis
          panel.spacing.y = unit(1.5, "lines"), #facet spacing for x axis
          strip.text.x = element_text(size=12, face="bold"), #facet labels
          strip.text.y = element_text(size=12, face="bold", angle = 270) #facet labels
    )
}

p_value_ext <- function(model,pt,cor) {
  if (pt == 1) {
    p <- glance(model)
    result <- paste("p = ",round(p$p.value,3))
  } else if (pt==2) {
    result <- paste("r = ",round(cor$estimate,3),"p = ",round(cor$p.value,3))
  }else if (pt==3) {
    result <- paste("ϱ = ",round(cor$estimate,3),"p = ",round(cor$p.value,3))
  } 
return(result)
}


