#packages for data vis and manipulating data
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
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


# stacked bar graphs
bar_stacked <- function(frac, yaxis, title)
{
  ggplot(Pfrac_NRP_stacked, aes(x=SiteName,y=Value, fill = Fraction)) +
    geom_bar(position = "stack", stat = "identity") +
    geom_errorbar(aes(ymax = Accumulation + Error, ymin = Accumulation - Error, width = 0.2)) + facet_grid(rows = vars(Area), cols = vars(SampleEvent)) + 
    ylab(yaxis) + scale_color_manual(values = rev(PNWColors::pnw_palette("Bay", length(table(frac$Fraction))))) + ggtitle(title)
}

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
