# R
#  plot average coverage from TSS file

args <-commandArgs(trailingOnly=TRUE)
input_file<-args[1]
output_file<-args[2]
text<-args[3]
ystart<-as.numeric(args[4])
ylimit<-as.numeric(args[5])

cov_data<-read.csv(input_file,header=TRUE,sep="\t")
png(filename =  output_file, width = 600, height = 600, units = "px", pointsize = 15, bg = "white", res = NA)
plot(cov_data$Position,cov_data$Mean.Cov,type="l",xlab="Position relative to TSS",ylab="Read depth",col="green",ylim=c(ystart,ylimit),bty="n")
polygon(c(cov_data$Position,rev(cov_data$Position)),c(cov_data$LowerBound,rev(cov_data$UpperBound)),col=rgb(215,25,28,51,maxColorValue=255), border = FALSE)
lines(cov_data$Position,cov_data$Mean.Cov,col=rgb(215,25,28,maxColorValue=255),lwd=3)
text(-500,0.2,text,col=rgb(215,25,28,maxColorValue=255))
dev.off()
