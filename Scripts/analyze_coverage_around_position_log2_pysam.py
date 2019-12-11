#! /usr/bin/env python

#Analyze read depth in comparison to transcription start

import sys
import argparse
import numpy
import scipy
import scipy.stats
import os.path
import pysam
import random


# Calculate mean and confidence intervals ###################################################################################
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*numpy.array(data)
    n = len(a)
    m, se = numpy.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

# Read in copy-number data ###################################################################################
def readCopyData(input_file):
    LOG2 = open(input_file,"r")
    header = LOG2.readline()
    log_list = list()
    for line in LOG2.readlines():
       chrom,start,end,log2 = line.rstrip().split("\t")
       log_list.append([chrom,int(start),int(end),float(log2)])
    LOG2.close()
    return log_list

# Calculate mean coverage in 10000bp upstream sequence ###################################################################################
def calcLocalMean(position, chrom, cnv_data):
    
    norm = None
    for region in cnv_data:
       reg_chrom,reg_start,reg_end,reg_log2 = region
       if chrom == reg_chrom and position > reg_start and position < reg_end:
           norm = numpy.power(2,reg_log2)
           break
    if norm:
        return norm
    else:
        return 1


# Calculate Coverage around position ###################################################################################
def calcCoverage(transcript_list,bam,cnv_data,args):
    sys.stderr.write("Coverage extraction started\n")
    sys.stderr.flush()
    position_lists=dict()

    #create a list of visited TSS not to count some more than once
    tss_visited = list()

    for i in range(-args.start,0):
        position_lists[i] = list()
    for i in range(0,args.end+1):
        position_lists[i] = list()

    line_count = 0
    skipped = 0
    used = 0
    for transcript in transcript_list:
        transcript_info = transcript.rstrip().split("\t")
        chrom = transcript_info[0]
        start = int(transcript_info[1])
        end = int(transcript_info[2])
        pos = (start+end)/2
        line_count += 1
        if (line_count  % 100 == 0):
            sys.stderr.write("\r "+str(line_count)+" regions analyzed. Skipped:"+str(skipped)+". Used: "+str(used))
            sys.stderr.flush()
        if len(transcript_info) < 6:
            forward = True
        elif transcript_info[5] == '+':
            forward = True
        else:
            forward = False
        if chrom+"_"+str(pos) in tss_visited:
            skipped += 1
            continue
        if pos-args.start < 1:
            skipped += 1
            continue
        if chrom.find("_") != -1:
            skipped += 1
            continue

        used += 1

        tss_visited.append(chrom+"_"+str(pos))
        coverage_tuple = bam.count_coverage(chrom, pos-args.start, pos+args.end+1, quality_threshold = args.mapq)
        coverage_list = list()
        for i in range(0,len(coverage_tuple[0])):
            coverage_list.append(coverage_tuple[0][i]+coverage_tuple[1][i]+coverage_tuple[2][i]+coverage_tuple[3][i])  
        #normlaize for copy-number variations
        if args.norm:
            normcov = calcLocalMean(pos, chrom,cnv_data) 
            if normcov == 0:
                normcov = 0.001       

        counter = 0
        for i in range(pos-args.start,pos+args.end+1):
            coverage = float(coverage_list[counter]) / args.mean_coverage
            if coverage > args.limit:
                continue
            if args.norm:
                coverage = coverage / normcov
            if forward:
                position_lists[i-pos].append(coverage)
            elif not forward:
                position_lists[-(i-pos)].append(coverage)
            counter += 1


    return(position_lists)
#######################################################################################################

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription start')
parser.add_argument('-bed','--bed-gene', dest='bed_file', 
                   help='BED file containing positions (chrom<TAB>start<TAB>end<TAB>[+,-])',required=True)
parser.add_argument('-b','--bam', dest='bam_file',
                   help='BAM file',required=True)
parser.add_argument('-s','--start', dest='start',
                   help='Start analyzing coverage at this point before region of interest [default:1000]',default=1000,type=int)
parser.add_argument('-e','--end', dest='end',
                   help='Stop analyzing coverage at this point after region of interest [default:1000]',default=1000,type=int)
parser.add_argument('-cov','--mean-coverage', dest='mean_coverage',
                   help='Mean coverage along the genome [default:1]',default=1,type=float)
parser.add_argument('-norm','--normalize', dest='norm',
                   help='Normalize by local coverage',action="store_true")
parser.add_argument('-norm-file','--normalize-file', dest='norm_log2',
                   help='Normalize by local copynumber from this file')
parser.add_argument('-limit','--coverage-limit', dest='limit',
                   help='discard coverage values above this threshold',default = 1000,type=float)
parser.add_argument('-m','--max-regions', dest='max_regions',
                   help='Use a maximum of this amount of regions',default=0,type=int)
parser.add_argument('-mapq','--mapping-quality-threshold', dest='mapq',
                   help='Only count coverage of reads with this mapping quality',default=0,type=int)


args = parser.parse_args()
sys.stderr.write("Bam file: "+args.bam_file+"\n")
sys.stderr.write("BED file: "+args.bed_file+"\n")
sys.stderr.write("Coverage: "+str(args.mean_coverage)+"\n")
sys.stderr.flush()

###############################################################################################
# Analyze data ###################################################################################

try:
    REFGENE = open(args.bed_file,"r")
except:
    print "Fail to open files specified"
    sys.exit(1)
    
#filter genes from genelist if specified
header = REFGENE.readline()
refgene_content = REFGENE.readlines()
target_content = refgene_content
if args.max_regions != 0 and args.max_regions < len(target_content):
    random.seed(args.bed_file)
    target_content=random.sample(target_content,args.max_regions)

sys.stderr.write("Starting analysis of "+str(len(target_content))+" regions"+"\n")
sys.stderr.flush()
#initialize input data
gene_count = 0
sys.stderr.write("\n")
sys.stderr.flush()

bam = pysam.AlignmentFile(args.bam_file, 'rb')

# read log2 data
cnv_data = readCopyData(args.norm_log2)

#collect all data
position_lists_all=calcCoverage(target_content,bam,cnv_data,args)

#output data
sys.stderr.write("--------------------------------------------------\n")
sys.stderr.write("\n"+str(len(position_lists_all[0]))+" TSS analyzed\n")
sys.stderr.flush()
print "Position\tMean Cov\tLowerBound\tUpperBound\tTSS analyzed\n"
for i in range(-args.start,args.end+1):
    if len(position_lists_all[i]) > 3:
        mean,lower_bound, upper_bound = mean_confidence_interval(position_lists_all[i])
        print str(i)+"\t"+str(mean)+"\t"+str(lower_bound)+"\t"+str(upper_bound)+"\t"+str(len(position_lists_all[i]))
    elif len(position_lists_all[i]) > 0:
        mean=numpy.mean(position_lists_all[i])
        print str(i)+"\t"+str(mean)+"\t"+str(len(position_lists_all[i]))
