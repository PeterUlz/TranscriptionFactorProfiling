#! R

args<-commandArgs(trailingOnly=TRUE)
infile<-args[1]
outfile<-args[2]

#analyze activity scores
library(ggplot2)
data<-read.csv(infile,header=FALSE,sep="\t")
#str(data)
#filtered_data<-subset(data, V3 < 100 & V3 > 0)
fit <- loess(V8 ~ log(V2),data=data)
data$prior<-predict(fit,data)

png(args[3])
ggplot()+geom_point(aes(x=data$V2,y=data$prior),color="blue")+geom_point(aes(x=data$V2,y=data$V8),color=rgb(0,0,0,0.3))+ylim(0,1)
dev.off()

data$adjusted_range<-data$V8 / data$prior
data$adjusted_range_correct_by_power <- data$adjusted_range * data$V5
data$adjusted_range_normalized <- data$adjusted_range / data$adjusted_range[which(data$V1 == "CTCF.50perc_hg19.tss")]
data$adjusted_range_rank <- rank(data$adjusted_range)
#ggplot(filtered_data,aes(x=V2,y=adjusted_range))+geom_point()

write.table(data,outfile,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

 
