#! R

library(WaveletComp)

findPeaks <- function (x, m = 3){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
       z <- i - m + 1
       z <- ifelse(z > 0, z, 1)
       w <- i + m + 1
       w <- ifelse(w < length(x), w, length(x))
       if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
     pks <- unlist(pks)
     pks
}

args<-commandArgs(trailingOnly=TRUE)

infile=args[1]
sample_name = args[2]
out_prefix = args[3]
all_out = args[4]

data<-read.csv(infile,header=TRUE,sep="\t")
n<-data$TSS.analyzed[1]
w<-analyze.wavelet(data,"Mean.Cov",n.sim=10)

png(paste0(out_prefix,".spectrum.png"))
wt.image(w,color.key="quantile")
dev.off()

png(paste0(out_prefix,".powers.png"))
period_powers = apply(w$Power,1,mean)
max_power = max(period_powers)
max_period = w$Period[which(period_powers==max_power)]
power_peaks = findPeaks(period_powers)
periods_of_peaks = w$Period[power_peaks]

power_peaks_in_range = power_peaks[which(periods_of_peaks > 135 & periods_of_peaks < 235)]
if (length(power_peaks_in_range) > 0) {
    max_power_in_range = max(period_powers[power_peaks_in_range])
    max_period_in_range = w$Period[which(period_powers==max_power_in_range)]
} else {
    max_power_in_range = 0
    max_period_in_range = 0
}
plot(w$Period,period_powers,pch=16,col="red",xlab="Period",ylab="Avg. power over range")
lines(w$Period,period_powers)
abline(v=max_period_in_range,col="red")
text(max_period_in_range+20,max_power_in_range,paste(format(round(max_period_in_range,2),nsmall=2),format(round(max_power_in_range,2),nsmall=2),sep=":"),adj=0)
if (max_period_in_range != max_period) {
    abline(v=max_period,col="blue")
    text(max_period+20,max_power,paste(format(round(max_period,2),nsmall=2),format(round(max_power,2),nsmall=2),sep=":"),adj=0)
}
dev.off()

png(paste0(out_prefix,".reconstructed.png"))
reconstr=reconstruct(w,sel.period=max_period_in_range)
dev.off()

range_reconstruct<-max(reconstr$series$Mean.Cov.r) - min(reconstr$series$Mean.Cov.r)
power_reconstruct<-sum((reconstr$series$Mean.Cov.r)**2)
sum_ampl_reconstruct<-sum(abs(reconstr$series$Mean.Cov.r))

write(paste(sample_name,n,max_period,max_power,max_period_in_range,max_power_in_range,range_reconstruct,power_reconstruct,sum_ampl_reconstruct,sep="\t"),file=all_out,append=TRUE)

