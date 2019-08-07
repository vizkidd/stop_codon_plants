#STOP specific mix model
library(plyr)

stops<-read.delim("Stop_codons",sep = " ")
stops <- stops[,colSums(is.na(stops))<nrow(stops)] ##CAN remove cols which are all NAs because they won't be in the tree and because we weren't able to download any sequence for the organism
#UAG=1
#UGA=2
#UAA=3
stop.proportions<-apply(stops,MARGIN = 1,function(x){
  UAG.count<-0
  UGA.count<-0
  UAA.count<-0
  for (ele in x) {
    switch (as.character(ele),
            "1" = {UAG.count=UAG.count+1},
            "2" = {UGA.count=UGA.count+1},
            "3" = {UAA.count=UAA.count+1}
    )
  }
  #stop.proportions<<-rbind(stop.proportions,
  return(c(UAG.count/length(x),UGA.count/length(x),UAA.count/length(x),(length(x)-(UAG.count+UGA.count+UAA.count))/length(x)))
})
rownames(stop.proportions)<-c("UAG","UGA","UAA","NA")
stop.proportions<-t(stop.proportions)
print("-----PROPORTIONS-----")
head(stop.proportions)

##IN TERMS of CLUSTERS
stop.counts<-apply(stops,MARGIN = 1,function(x){
  UAG.count<-0
  UGA.count<-0
  UAA.count<-0
  for (ele in x) {
    switch (as.character(ele),
            "1" = {UAG.count=UAG.count+1},
            "2" = {UGA.count=UGA.count+1},
            "3" = {UAA.count=UAA.count+1}
    )
  }
  #stop.proportions<<-rbind(stop.proportions,
  return(c(UAG.count,UGA.count,UAA.count,(length(x)-(UAG.count+UGA.count+UAA.count))))
})
rownames(stop.counts)<-c("UAG","UGA","UAA","NA")
stop.counts<-t(stop.counts)
print("-----COUNTS-----")
head(stop.counts)

org.stop.freq<-t(apply(stops, MARGIN =c(2), FUN=function(x){cnt<-count(x); cnt<-cnt[,-1]; cnt}))
head(org.stop.freq)
colnames(org.stop.freq)<-c("UAG","UGA","UAA","NA")
org.stop.freq<-as.data.frame(org.stop.freq)
colSums(org.stop.freq)
colMeans(org.stop.freq)

##choose BEST represented org
chosen.org<-rownames(org.stop.freq)[which(org.stop.freq[,4]==min(org.stop.freq[,4]))]
print("-----BEST REPRESENTED ORGANISM-----")
chosen.org ##BEST represented org
chosen.org.stops<- stops[,chosen.org]
#CONVERT NAs to '-'
chosen.org.stops[is.na(chosen.org.stops)]<- "-"
chosen.clust.UAG<-rownames(stops)[which(chosen.org.stops=="1")]
chosen.clust.UGA<-rownames(stops)[which(chosen.org.stops=="2")]
chosen.clust.UAA<-rownames(stops)[which(chosen.org.stops=="3")]

mparams<-read.table("Model_parameters.mgf1x4")
trees<-read.delim("Trees.orig",header = FALSE)
colnames(trees)<- c("clust","tree")
extMG<-read.delim("extMG_params",sep=" ",header = FALSE)

##WRITE files into seperate dirs to apply mix model on seperate stops

##SAVE STOPS of represented organism in which_stop.txt
which.stops<-data.frame(Stops=stops[,chosen.org],row.names = rownames(stops))
#which.stops<- na.omit(which.stops)
which.stops[which(is.na(which.stops[,1])),]<- "All"
which.stops[which.stops[,1]==1,]<- "UAG"
which.stops[which.stops[,1]==2,]<- "UGA"
which.stops[which.stops[,1]==3,]<- "UAA"
write.table(which.stops,"which_stop.txt",sep = " ",quote = FALSE,col.names = FALSE)

##UAG

extMG.UAG<- extMG[which(extMG$V9 %in% chosen.clust.UAG),]
stops.UAG<- stops[which(rownames(stop.counts) %in% chosen.clust.UAG),]
trees.UAG<- trees[which(trees$clust %in% chosen.clust.UAG),]
mparams.UAG<- mparams[which(mparams$V1 %in% chosen.clust.UAG),]
write.table(trees.UAG$tree,"UAG/Trees",sep = " ",quote= FALSE,col.names = FALSE,row.names = FALSE)
write.table(mparams.UAG,"UAG/Model_parameters.mgf1x4",sep = " ",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(stops.UAG,"UAG/Stop_codons",quote = FALSE)
write.table(extMG.UAG,"UAG/extMG",quote = FALSE)

##UGA

extMG.UGA<- extMG[which(extMG$V9 %in% chosen.clust.UGA),]
stops.UGA<- stops[which(rownames(stop.counts) %in% chosen.clust.UGA),]
trees.UGA<- trees[which(trees$clust %in% chosen.clust.UGA),]
mparams.UGA<- mparams[which(mparams$V1 %in% chosen.clust.UGA),]
write.table(trees.UGA$tree,"UGA/Trees",sep = " ",quote= FALSE,col.names = FALSE,row.names = FALSE)
write.table(mparams.UGA,"UGA/Model_parameters.mgf1x4",sep = " ",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(stops.UGA,"UGA/Stop_codons",quote = FALSE)
write.table(extMG.UGA,"UGA/extMG",quote = FALSE)

##UAA

extMG.UAA<- extMG[which(extMG$V9 %in% chosen.clust.UAA),]
stops.UAA<- stops[which(rownames(stop.counts) %in% chosen.clust.UAA),]
trees.UAA<- trees[which(trees$clust %in% chosen.clust.UAA),]
mparams.UAA<- mparams[which(mparams$V1 %in% chosen.clust.UAA),]
write.table(trees.UAA$tree,"UAA/Trees",sep = " ",quote= FALSE,col.names = FALSE,row.names = FALSE)
write.table(mparams.UAA,"UAA/Model_parameters.mgf1x4",sep = " ",quote = FALSE,col.names = FALSE,row.names = FALSE)
write.table(stops.UAA,"UAA/Stop_codons",quote = FALSE)
write.table(extMG.UAA,"UAA/extMG",quote = FALSE)
