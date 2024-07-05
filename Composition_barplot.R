library(ggplot2)
library(xlsx)
library(forcats)

###Split TSS distance to groups and plot

Dists <- read.csv('./CHIP/Final_ann_V_WT.csv')
Dists <- Dists[Dists$Fold<(-log2(1.5)),]
names(Dists)[1]<-"PeakID"
Dists$AbsDistTSS <- abs(Dists$distanceToTSS)
Dists$DistGroup <- cut(Dists$AbsDistTSS,
                                  breaks = c(-1,100,1000,10000,100000,Inf),
                                  labels = c(c("0-100","100-1000","1000-10000","10000-100000",">100000")))

Dists <- data.frame(Dists$PeakID,Dists$DistGroup)
names(Dists)<-c("PeakID","DistGroup")
Dists$DistGroup2 <- fct_relevel(Dists$DistGroup,c(">100000","10000-100000","1000-10000",
                             "100-1000","0-100"))
Dists$type <- "V2084D-WT Down"

Dists2 <- read.csv('./CHIP/Res_Vsh183_allpeaks.csv')
names(Dists2)[1]<-"PeakID"
Dists2$AbsDistTSS <- abs(Dists2$distanceToTSS)
Dists2$DistGroup <- cut(Dists2$AbsDistTSS,
                         breaks = c(-1,100,1000,10000,100000,Inf),
                         labels = c(c("0-100","100-1000","1000-10000","10000-100000",">100000")))

Dists2 <- data.frame(Dists2$PeakID,Dists2$DistGroup)
names(Dists2)<-c("PeakID","DistGroup")
Dists2$DistGroup2 <- fct_relevel(Dists2$DistGroup,c(">100000","10000-100000","1000-10000",
                                                  "100-1000","0-100"))
Dists2$type <- "VshSTUB1-V2084D Rescued"
DistsFinal <- rbind(Dists,Dists2)

##Percentage stacked bar
p<-ggplot(DistsFinal, aes(x=type,fill = DistGroup2)) +
  geom_bar(width = 0.5, position = "fill") + 
  scale_fill_brewer(palette = "Blues",direction = -1) +   
  scale_y_continuous(labels = scales::percent) +  
  guides(fill=guide_legend(title = "Distance To TSS")) +
  labs(x = " ",
       y = "peaks") +
  theme(axis.text.x = element_text (angle = 45, 
                                    vjust = 1, 
                                    hjust=1),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(),
        panel.grid.major.y = element_line(colour ="#EEEEEE"))
p

