#! /usr/bin/python
import subprocess
import os
import operator
###########################################################################################
# Split ChIP-Seq reference data into tf Class
def convert(factor):
    infilename = factor.split(" ")[0].replace("/","-")+".hg38.bed"
    outfilename = factor.split(" ")[0].replace("/","-")+".hg19.bed"
    subprocess.call(["liftOver","-bedPlus=6","-tab",infilename,"hg38ToHg19.over.chain",outfilename,"unmapped"])
    subprocess.call(["rm",infilename])
############################################################################################
def filterFactor(factor):
    infilename = factor.split(" ")[0].replace("/","-")+".hg19.bed"
    outfilename = factor.split(" ")[0].replace("/","-")+".Top1000sites.bed"
    IN = open(infilename,"r")
    file_dict = dict()
    lines = IN.readlines()
    for line in lines:
        chrom,start,end,tf,peaks,strand = line.rstrip().split("\t")
        file_dict[line] = int(peaks)
    sorted_file = sorted(file_dict.items(), key=operator.itemgetter(1))
    sites = len(sorted_file)
    if sites < 1000:
        print "Too few sites"
        return
    OUT = open(outfilename,"w")
    for i in range(1,1001):
        OUT.write(sorted_file[sites-i][0])
    OUT.close()
    IN.close()
############################################################################################  
IN = open("../GTRD/human_meta_clusters.interval","r")
line = IN.readline()
this_factor = None
OUT = open("stats","w")
while line:
    if line[0] == "#":
        line = IN.readline()
        continue
    chrom,start,end,summit,tfclass,name,cells,treatment,exp,peak,peak_count,exp_count,peak_count = line.rstrip().split("\t")
    if chrom == "chrX" or chrom == "chrY" or "_" in chrom:
        line = IN.readline()
        continue
    position = int(start)+int(summit)
    new_end = position+1
    if name != this_factor:
        OUT.close()
        if this_factor != None:
            convert(this_factor)
            filterFactor(this_factor)
        this_factor = name
        OUT = open(name.split(" ")[0].replace("/","-")+".hg38.bed","w")
        OUT.write('\t'.join([chrom,str(position),str(new_end),name,peak_count,"+"])+"\n")
    else:
        OUT.write('\t'.join([chrom,str(position),str(new_end),name,peak_count,"+"])+"\n")
    line = IN.readline()
OUT.close()




