# R
#  plot average coverage from TSS file

args <-commandArgs(trailingOnly=TRUE)
input_file1<-args[1]
input_file2<-args[2]
input_file3<-args[3]
output_file<-args[4]
ystart<-as.numeric(args[5])
ylimit<-as.numeric(args[6])

cov_data1<-read.csv(input_file1,header=TRUE,sep="\t")
cov_data2<-read.csv(input_file2,header=TRUE,sep="\t")
cov_data3<-read.csv(input_file3,header=TRUE,sep="\t")
png(filename =  output_file, width = 600, height = 600, units = "px", pointsize = 15, bg = "white", res = NA)
plot(cov_data1$Position,cov_data1$Mean.Cov,type="l",xlab="Position relative to TSS",ylab="Read depth",col="green",ylim=c(ystart,ylimit),bty="n")
polygon(c(cov_data1$Position,rev(cov_data1$Position)),c(cov_data1$LowerBound,rev(cov_data1$UpperBound)),col=rgb(215,25,28,51,maxColorValue=255), border = FALSE)
lines(cov_data1$Position,cov_data1$Mean.Cov,col=rgb(215,25,28,maxColorValue=255),lwd=3)
polygon(c(cov_data2$Position,rev(cov_data2$Position)),c(cov_data2$LowerBound,rev(cov_data2$UpperBound)),col=rgb(44,123,182,51,maxColorValue=255), border = FALSE)
lines(cov_data2$Position,cov_data2$Mean.Cov,col=rgb(44,123,182,maxColorValue=255),lwd=3)
polygon(c(cov_data3$Position,rev(cov_data3$Position)),c(cov_data3$LowerBound,rev(cov_data3$UpperBound)),col=rgb(44,182,30,51,maxColorValue=255), border = FALSE)
lines(cov_data3$Position,cov_data3$Mean.Cov,col=rgb(44,182,30,maxColorValue=255),lwd=3)
dev.off()
