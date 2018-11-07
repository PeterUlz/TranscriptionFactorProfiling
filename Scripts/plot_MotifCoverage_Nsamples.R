# R
#  plot average coverage from TSS file


args <-commandArgs(trailingOnly=TRUE)
input_dir<-args[1]
ystart<-as.numeric(args[2])
ylimit<-as.numeric(args[3])
output_file<-args[4]
files=Sys.glob(paste0(input_dir,"*.tss"))
colors=rainbow(length(files))
cov_data1<-read.csv(files[1],header=TRUE,sep="\t")
png(filename =  output_file, width = 1200, height = 600)
plot(cov_data1$Position,cov_data1$Mean.Cov,type="n",xlab="Position relative to TSS",ylab="Read depth",col="green",ylim=c(ystart,ylimit),bty="n")

for (i in 1:length(files)) {
    cov_data1<-read.csv(files[i],header=TRUE,sep="\t")
    lines(cov_data1$Position,cov_data1$Mean.Cov,col=colors[i],ylim=c(ystart,ylimit),bty="n")
    polygon(c(cov_data1$Position,rev(cov_data1$Position)),c(cov_data1$LowerBound,rev(cov_data1$UpperBound)),col=rgb(t(col2rgb(colors[i])),alpha=51,maxColorValue=255), border = FALSE)
}
legend("bottom",files,ncol=5,fill=colors)
dev.off()

