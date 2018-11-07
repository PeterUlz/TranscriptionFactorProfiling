#! R

args<-commandArgs(trailingOnly=TRUE)
infile<-args[1]
outfile<-args[2]

#analyze activity scores
library(ggplot2)
data<-read.csv(infile,header=FALSE,sep="\t")

data$adjusted_range_rank <- rank(data$V3)
#ggplot(filtered_data,aes(x=V2,y=adjusted_range))+geom_point()

write.table(data,outfile,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

 
