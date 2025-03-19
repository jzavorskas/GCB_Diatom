# Setup
library(ggplot2)
library(greekLetters)
library(gridExtra)
setwd("C:/Users/josep/miniconda3/envs/dFBA/!!! Notebooks")
Fig1Data <- read.csv("./surfplotallMayTranspose.csv",
                     header = FALSE)
Fig2Data <- read.csv("./SuccessionMarchTranspose.csv",
                     header = FALSE)
Fig6aData <- read.csv("./Fig6AData.csv",
                      header = FALSE)
Fig6bData <- read.csv("./Fig6BData.csv",
                      header = FALSE)
Fig6cData <- read.csv("./Fig6CData.csv",
                      header = FALSE)
Fig6dData <- read.csv("./Fig6DData.csv",
                      header = FALSE)

### Figure 1

ggplot(Fig1Data, aes(x=V1, y=V3)) + theme_bw() +
       geom_point(aes(col=as.factor(V2)), size=0.75) +
       scale_shape_manual(values=c(55, 65, 75, 85)) + 
       scale_color_manual(values=c('#fcaa53','#cc6a04','#56B4E9','#000080')) +
       labs(y=expression(paste("Thalassiosira ", italic("sp. "), "Biomass (mg dry weight/L)")),
            x="Time (hours)",color="degrees\nNorth\nLatitude") +
       theme(axis.text.x = element_text(face="bold",size=9),
            axis.text.y = element_text(face="bold",size=9),
            axis.title = element_text(size=10.5),
            legend.text = element_text(face="bold",size=10),
            legend.title = element_text(face="bold",size=9,hjust=0.5),
            legend.background = element_rect(fill = "#fbfbfb",colour='black'),
            legend.key = element_rect(fill = "#fbfbfb")) +
        guides(color = guide_legend(override.aes = list(size=3)))
             
ggsave(filename = "Fig1GG.png", width=5, height=3.5, device='png', dpi=300)                                     
   
### Figure 2a

Fig2a <- ggplot(Fig1Data, aes(x=V1, y=V5)) + theme_bw() +
          geom_point(aes(col=as.factor(V2)), size=0.75) +
          scale_shape_manual(values=c(55, 65, 75, 85)) + 
          scale_color_manual(values=c('#fcaa53','#cc6a04','#56B4E9','#000080')) +
          ylab(expression("NO"[3] ~ Concentration ~ paste("\n (", mu, "mol/L)"))) +
          xlab("Time (hours)") +
          labs(tag="A") +
          #labs(color="degrees\nNorth\nLatitude",tag="A") +
          theme(axis.text.x = element_text(face="bold",size=9),
                axis.text.y = element_text(face="bold",size=9),
                axis.title = element_text(size=10.5),
                legend.position = "none") +
                #legend.text = element_text(face="bold",size=10),
                #legend.title = element_text(face="bold",size=9,hjust=0.5),
                #legend.background = element_rect(fill = "#fbfbfb",colour='black'),
                #legend.key = element_rect(fill = "#fbfbfb")) +
          # guides(color = guide_legend(override.aes = list(size=3))) +
          scale_y_continuous(limits = c(0, 10),breaks=c(0,3,6,9))

# ggsave(filename = "Fig2aGG.png", width=5, height=3.5, device='png', dpi=300) 

## Figure 2b       

Fig2b <- ggplot(Fig1Data, aes(x=V1, y=V6)) + theme_bw() +
          geom_point(aes(col=as.factor(V2)), size=0.75) +
          scale_shape_manual(values=c(55, 65, 75, 85)) + 
          scale_color_manual(values=c('#fcaa53','#cc6a04','#56B4E9','#000080')) +
          ylab(expression("Si(OH)"[4] ~ Concentration ~ paste("\n (", mu, "mol/L)"))) +
          xlab("Time (hours)") +
          labs(color="degrees\nNorth\nLatitude",tag="B") +
          theme(axis.text.x = element_text(face="bold",size=9),
                axis.text.y = element_text(face="bold",size=9),
                axis.title = element_text(size=10.5),
                legend.text = element_text(face="bold",size=10),
                legend.title = element_text(face="bold",size=9,hjust=0.5),
                legend.background = element_rect(fill = "#fbfbfb",colour='black'),
                legend.key = element_rect(fill = "#fbfbfb")) +
          guides(color = guide_legend(override.aes = list(size=3))) +
          scale_y_continuous(limits = c(1, 4),breaks=c(0,1,2,3,4),minor_breaks=NULL) +
          geom_hline(yintercept=2, linetype="21", color = "black",linewidth=1) +
          annotate("text", x=1515, y=2.125, 
                   label = expression(Expected ~ Final ~ "Si(OH)"[4] ~ Concentration),
                   size = unit(3, "pt"))

Fig2Combo <- grid.arrange(Fig2a, Fig2b, nrow = 1, widths=c(0.9,1))

ggsave(filename = "Fig2SingleLegendGG.png", plot=Fig2Combo, width=9, height=3.5, device='png', dpi=300)

## Figure 4

ggplot(Fig2Data, aes(x=V1, y=V4)) + theme_bw() +
  geom_point(aes(col=as.factor(V2)), size=0.75) +
  scale_shape_manual(values=c(55, 65, 75, 85)) + 
  scale_color_manual(values=c('#fcaa53','#cc6a04','#56B4E9','#000080')) +
  labs(y=expression(paste("Chaetoceros ", italic("sp. "), "Biomass (mg dry weight/L)")),
       x="Time (hours)",color="degrees\nNorth") +
  theme(axis.text.x = element_text(face="bold",size=9),
        axis.text.y = element_text(face="bold",size=9),
        axis.title = element_text(size=10.5),
        legend.text = element_text(face="bold",size=10),
        legend.title = element_text(face="bold",size=9,hjust=0.5),
        legend.background = element_rect(fill = "#fbfbfb",colour='black'),
        legend.key = element_rect(fill = "#fbfbfb")) +
  guides(color = guide_legend(override.aes = list(size=3)))

ggsave(filename = "Fig4GG.png", width=5, height=3.5, device='png', dpi=300)

## Figure 5

ggplot(Fig2Data, aes(x=V1, y=(V3+V4))) + theme_bw() + 
  geom_point(aes(col=as.factor(V2)), size=0.75) +
  scale_shape_manual(values=c(55, 65, 75, 85)) + 
  scale_color_manual(values=c('#fcaa53','#cc6a04','#56B4E9','#000080')) +
  labs(y=expression(paste("Combined Diatom Biomass (mg dry weight/L)")),
       x="Time (hours)",color="degrees\nNorth") +
  theme(axis.text.x = element_text(face="bold",size=9),
        axis.text.y = element_text(face="bold",size=9),
        axis.title = element_text(face="bold",size=9),
        legend.text = element_text(face="bold",size=10),
        legend.title = element_text(face="bold",size=9,hjust=0.5),
        legend.background = element_rect(fill = "#fbfbfb",colour='black'),
        legend.key = element_rect(fill = "#fbfbfb")) +
  guides(color = guide_legend(override.aes = list(size=2.5)))

ggsave(filename = "Fig5GG.png", width=5, height=3.5, device='png', dpi=300)

## Figure 6
# Four subplots, one for each latitude
# Colors by nutrient type, not latitude
# Need to make four .csv files, one for each lat (easiest way)

cols <- as.list(setNames(c('#fcaa53','#A14700','#418CD9'),
                         c("Nitrate","Silicon","Iron"))) 
# labls <- list(expression("NO"[3]), expression("Si(OH)"[4]),"dFe")
bks <- c("Nitrate","Silicon","Iron")

Fig6a <- ggplot(Fig6aData, aes(x=V1)) + theme_bw() +
  geom_point(aes(x=V1,y=V8,col=bks[1]), size=0.75) +
  geom_point(aes(x=V1,y=V9,col=bks[2]), size=0.75) +
  geom_point(aes(x=V1,y=V10/0.001,col=bks[3]), size=0.75) +
  scale_color_manual(values=cols,breaks=bks) +
  ylab(expression("NO"[3] ~ or ~ "Si(OH)"[4] ~ Concentration ~
                  paste("\n (", mu, "mol/L)"))) +
  xlab("Time (hours)") +
  labs(color="Nutrient",tag="A") +
  scale_y_continuous(limits = c(0, 10),breaks=c(0,2,4,6,8,10),minor_breaks=NULL,
                     sec.axis = sec_axis(transform=~.*1, 
                        name=paste("dFe Concentration (nmol/L)"),
                        breaks=c(0,2,4,6,8,10))) +
  theme(axis.text.x = element_text(face="bold",size=9),
        axis.text.y = element_text(face="bold",size=9),
        axis.text.y.right = element_text(color = '#418CD9'),
        axis.title = element_text(size=10.5),
        axis.title.y.right = element_text(color = '#418CD9'),
        legend.text = element_text(size=10),
        legend.title = element_text(face="bold",size=9,hjust=0.5),
        legend.background = element_rect(fill = "#fbfbfb",colour='black'),
        legend.key = element_rect(fill = "#fbfbfb"),
        legend.position = "inside",
        legend.position.inside=c(0.85,0.775),) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  geom_hline(yintercept=2, linetype="21", color = "black",linewidth=1) +
  annotate("text", x=3100, y=3, 
           label = expression("Expected Final"),
           size = unit(3, "pt"),colour="#A14700") +
  annotate("text", x=3100, y=2.4, 
          label = expression("\nSi(OH)"[4] ~ Concentration),
           size = unit(3, "pt"),colour="#A14700")

Fig6b <- ggplot(Fig6bData, aes(x=V1)) + theme_bw() +
  geom_point(aes(x=V1,y=V8,col=bks[1]), size=0.75) +
  geom_point(aes(x=V1,y=V9,col=bks[2]), size=0.75) +
  geom_point(aes(x=V1,y=V10/0.001,col=bks[3]), size=0.75) +
  scale_color_manual(values=cols,breaks=bks) +
  ylab(expression("NO"[3] ~ or ~ "Si(OH)"[4] ~ Concentration ~
                    paste("\n (", mu, "mol/L)"))) +
  xlab("Time (hours)") +
  labs(color="Nutrient",tag="B") +
  scale_y_continuous(limits = c(0, 10),breaks=c(0,2,4,6,8,10),minor_breaks=NULL,
                     sec.axis = sec_axis(transform=~.*1, 
                                         name=paste("dFe Concentration (nmol/L)"),
                                         breaks=c(0,2,4,6,8,10))) +
  theme(axis.text.x = element_text(face="bold",size=9),
        axis.text.y = element_text(face="bold",size=9),
        axis.text.y.right = element_text(color = '#418CD9'),
        axis.title = element_text(size=10.5),
        axis.title.y.right = element_text(color = '#418CD9'),
        legend.text = element_text(size=10),
        legend.title = element_text(face="bold",size=9,hjust=0.5),
        legend.background = element_rect(fill = "#fbfbfb",colour='black'),
        legend.key = element_rect(fill = "#fbfbfb"),
        legend.position = "inside",
        legend.position.inside=c(0.85,0.775)) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  geom_hline(yintercept=2, linetype="21", color = "black",linewidth=1) +
  annotate("text", x=3100, y=3, 
           label = expression("Expected Final"),
           size = unit(3, "pt"),colour="#A14700") +
  annotate("text", x=3100, y=2.4, 
           label = expression("\nSi(OH)"[4] ~ Concentration),
           size = unit(3, "pt"),colour="#A14700")

Fig6c <- ggplot(Fig6cData, aes(x=V1)) + theme_bw() +
  geom_point(aes(x=V1,y=V8,col=bks[1]), size=0.75) +
  geom_point(aes(x=V1,y=V9,col=bks[2]), size=0.75) +
  geom_point(aes(x=V1,y=V10/0.001,col=bks[3]), size=0.75) +
  scale_color_manual(values=cols,breaks=bks) +
  ylab(expression("NO"[3] ~ or ~ "Si(OH)"[4] ~ Concentration ~
                    paste("\n (", mu, "mol/L)"))) +
  xlab("Time (hours)") +
  labs(color="Nutrient",tag="C") +
  scale_y_continuous(limits = c(0, 10),breaks=c(0,2,4,6,8,10),minor_breaks=NULL,
                     sec.axis = sec_axis(transform=~.*1, 
                                         name=paste("dFe Concentration (nmol/L)"),
                                         breaks=c(0,2,4,6,8,10))) +
  theme(axis.text.x = element_text(face="bold",size=9),
        axis.text.y = element_text(face="bold",size=9),
        axis.text.y.right = element_text(color = '#418CD9'),
        axis.title = element_text(size=10.5),
        axis.title.y.right = element_text(color = '#418CD9'),
        legend.text = element_text(size=10),
        legend.title = element_text(face="bold",size=9,hjust=0.5),
        legend.background = element_rect(fill = "#fbfbfb",colour='black'),
        legend.key = element_rect(fill = "#fbfbfb"),
        legend.position = "inside",
        legend.position.inside=c(0.85,0.775)) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  geom_hline(yintercept=2, linetype="21", color = "black",linewidth=1) +
  annotate("text", x=3100, y=3, 
           label = expression("Expected Final"),
           size = unit(3, "pt"),colour="#A14700") +
  annotate("text", x=3100, y=2.4, 
           label = expression("\nSi(OH)"[4] ~ Concentration),
           size = unit(3, "pt"),colour="#A14700")

Fig6d <- ggplot(Fig6dData, aes(x=V1)) + theme_bw() +
  geom_point(aes(x=V1,y=V8,col=bks[1]), size=0.75) +
  geom_point(aes(x=V1,y=V9,col=bks[2]), size=0.75) +
  geom_point(aes(x=V1,y=V10/0.001,col=bks[3]), size=0.75) +
  scale_color_manual(values=cols,breaks=bks) +
  ylab(expression("NO"[3] ~ or ~ "Si(OH)"[4] ~ Concentration ~
                    paste("\n (", mu, "mol/L)"))) +
  xlab("Time (hours)") +
  labs(color="Nutrient",tag="D") +
  scale_y_continuous(limits = c(0, 10),breaks=c(0,2,4,6,8,10),minor_breaks=NULL,
                     sec.axis = sec_axis(transform=~.*1, 
                                         name=paste("dFe Concentration (nmol/L)"),
                                         breaks=c(0,2,4,6,8,10))) +
  theme(axis.text.x = element_text(face="bold",size=9),
        axis.text.y = element_text(face="bold",size=9),
        axis.text.y.right = element_text(color = '#418CD9'),
        axis.title = element_text(size=10.5),
        axis.title.y.right = element_text(color = '#418CD9'),
        legend.text = element_text(size=10),
        legend.title = element_text(face="bold",size=9,hjust=0.5),
        legend.background = element_rect(fill = "#fbfbfb",colour='black'),
        legend.key = element_rect(fill = "#fbfbfb"),
        legend.position = "inside",
        legend.position.inside=c(0.85,0.775)) +
  guides(color = guide_legend(override.aes = list(size=3))) +
  geom_hline(yintercept=2, linetype="21", color = "black",linewidth=1) +
  annotate("text", x=3100, y=3, 
           label = expression("Expected Final"),
           size = unit(3, "pt"),colour="#A14700") +
  annotate("text", x=3100, y=2.4, 
           label = expression("\nSi(OH)"[4] ~ Concentration),
           size = unit(3, "pt"),colour="#A14700")

Fig6Combo <- grid.arrange(Fig6a, Fig6b, Fig6c, Fig6d, nrow = 2, ncol = 2)

ggsave(plot=Fig6Combo, filename = "Fig6GG.png", width=10, height=7, device='png', dpi=300)