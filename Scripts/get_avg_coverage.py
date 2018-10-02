#! /usr/bin/python

import sys

chromosomes = dict()
chrom_sizes = dict()
IN = open(sys.argv[1],"r")

for line in IN.readlines():
    chrom,cov,bases,size,fraction = line.rstrip().split("\t")
    if chrom not in chromosomes.keys():
        chromosomes[chrom] = 0
        chrom_sizes[chrom] = int(size)
    chromosomes[chrom] += int(cov)*int(bases)

autosome_coverage = 0
auosome_size = 0
for autosome in range(1,23):
    chromosome = "chr"+str(autosome)
    autosome_size += chrom_sizes[chromosome]
    autosome_coverage += chromosomes[chromosome] * chrom_sizes[chromosome]
autosome_coverage = float(autosome_coverage) / float(autosome_size)

for chrom in chromosomes.keys():
    print chrom+"\t"+str(float(chromosomes[chrom])/float(chrom_sizes[chrom]))
    print "Autosomes\t"+str()
