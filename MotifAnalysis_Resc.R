library(limma)
library(ggplot2)

#Plot logP

motifs <- read.delim('../Processed/HomerMotif/Res_Vsh_V/knownResults.txt')
motifs$Name <- strsplit2(motifs$Motif.Name,'\\/')[,1]
motifs$logP <- (-motifs$Log.P.value)
motifs <- motifs[order(motifs$Log.P.value),]
motifs$q.value <- motifs$q.value..Benjamini.
motifstop <- motifs[1:10,]


p <- ggplot(motifstop, aes(x=reorder(Name,logP),y=logP,fill="red"))+
  geom_bar(stat = "identity")+
  coord_flip()+
  xlab("Motif")+
  ylab("-log(p)")+
  ggtitle("Rescued peaks")+
  #scale_fill_gradient(low="red",high = "pink")+
  theme_classic()

p

