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

for chrom in chromosomes.keys():
    print chrom+"\t"+str(float(chromosomes[chrom])/float(chrom_sizes[chrom]))

