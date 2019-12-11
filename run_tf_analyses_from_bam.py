#! /usr/bin/env python

# Analyze all possible things from BAM-file

import sys
import argparse
from subprocess import call
import numpy
import scipy
import scipy.stats
import os.path
import os
import glob

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze epigenetic traces in cfDNA')
parser.add_argument('-b','--bam', dest='bam_file',
                   help='BAM file',required=True)
parser.add_argument('-o','--output', dest='name',
                   help='Output name for files and directory',required=True)
parser.add_argument('-cov','--mean-coverage', dest='mean_coverage',
                   help='Mean coverage along the genome [default:1]',default=1,type=float)
parser.add_argument('-ylimit','--plot-y-limit', dest='ylimit',
                   help='Plotting until this limit on y-axis [default:1.5]',default=1.5,type=float)
parser.add_argument('-norm-file','--normalize-file', dest='norm_log2',
                   help='Normalize by local copynumber from this file')
parser.add_argument('-calccov','--calculate-mean-coverage', dest='calc_cov',
                   help='Specify whether genome read depths should be calculated',action="store_true")
parser.add_argument('-hg38','--hg38', dest='hg38',
                   help='Use hg38 coordinates [default: hg19]',action="store_true")
parser.add_argument('-a','--analysis', dest='analysis',
                   help='Specify type of analysis (all|enhancer|histone|tf|ctcf|...)',required=True)
parser.add_argument('-tf','--trans-factor', dest='tf',
                   help='Specify transcription factor for VirChip data')

args = parser.parse_args()
####################################################################################################
# setup structure
print "Setup structure"
if not os.path.isdir(args.name):
    os.mkdir(args.name)

####################################################################################################
# get genomewide coverage from bedtools genomecoverage
if args.calc_cov:
    print "Calc avg. coverage"
    OUTPUT=open(args.name.rstrip("/")+"/"+args.name+".coverage","w")
    if args.hg38:
        call(["./Software/bedtools","genomecov","-ibam",args.bam_file,"-g","./Ref/hg38.chrom_sizes.txt"],stdout=OUTPUT)
    else:
        call(["./Software/bedtools","genomecov","-ibam",args.bam_file,"-g","./Ref/hg19.chrom_sizes.txt"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT=open(args.name.rstrip("/")+"/"+args.name+".short_coverage","w")
    call(["./Scripts/get_avg_coverage.py",args.name.rstrip("/")+"/"+args.name+".coverage"],stdout=OUTPUT)
    OUTPUT.close()
    INPUT = open(args.name.rstrip("/")+"/"+args.name+".short_coverage","r")
    avg_coverage = 1
    for line in INPUT.readlines():
        chrom,cov = line.rstrip().split("\t")
        if chrom == "genome":
            avg_coverage = cov
    INPUT.close()
else:
    print "Skipping genomewide-coverage calculation using mean coverage: "+str(args.mean_coverage)
    avg_coverage = args.mean_coverage
####################################################################################################
# print statistics:
print "Write Logs"
OUT=open(args.name.rstrip("/")+"/log.txt","w")
OUT.write("BAM:\t"+args.bam_file+"\n")
OUT.write("Norm File:\t"+args.norm_log2+"\n")
OUT.write("cov:\t"+str(avg_coverage)+"\n")
OUT.write("analysis:\t"+args.analysis+"\n")
OUT.close()
####################################################################################################
# get chromosome coverage from output of bedtools genomecoverage
def getChromCoverage(chromosome,args):
    print args.name.rstrip("/")+"/"+args.name+".short_coverage"
    if not os.path.isfile(args.name.rstrip("/")+"/"+args.name+".short_coverage"):
        print "Coverage file not found"
        sys.exit(1)
    INPUT = open(args.name.rstrip("/")+"/"+args.name+".short_coverage","r")
    avg_coverage = 1
    found = False
    for line in INPUT.readlines():
        chrom,cov = line.rstrip().split("\t")
        if chrom == chromosome:
            avg_coverage = cov
            found = True
    INPUT.close()
    if found:
        return avg_coverage
    else:
        print "Chromosome not found"
        sys.exit(1)
####################################################################################################
# CTCF analysis
def ctcf(args,avg_coverage):
    print "Analyze CTCF sites"
    if not os.path.isdir(args.name.rstrip("/")+"/CTCF"):
        os.mkdir(args.name.rstrip("/")+"/CTCF")
    #OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_In_Insulated_Neighbourhoods.tss","w")
    #call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/FIMO_ChIP_CTCF_at_Insulated_Neighbourhoods.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    #OUTPUT.close()
    #OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_Outside_Insulated_Neighbourhoods.tss","w")
    #call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/FIMO_ChIP_CTCF_outside_Insulated_Neighbourhoods.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    #OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_GTRD_In_Insulated_Neighbourhoods.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/CTCF.GTRD.Insulated.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-m","100000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_GTRD_Outside_Insulated_Neighbourhoods.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/CTCF.GTRD.NonInsulated.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-m","100000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_GTRD_50perc_In_Insulated_Neighbourhoods.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/CTCF.GTRD.50perc.Insulated.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-m","100000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_GTRD_50perc_Outside_Insulated_Neighbourhoods.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/CTCF.GTRD.50perc.NonInsulated.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-m","100000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_ultraconserved.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/Ultraconserved_CTCF.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-m","100000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_GTRD_proximalTSS.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/CTCF.hg19.sorted.bed.proximal.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-m","100000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_GTRD_distalTSS.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/CTCF.hg19.sorted.bed.distal.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-m","100000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_GTRD_50perc_proximalTSS.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/CTCF.50perc_hg19.sorted.bed.proximal.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-m","100000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/CTCF"+"/CTCF_GTRD_50perc_distalTSS.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/CTCF/CTCF.50perc_hg19.sorted.bed.distal.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-m","100000"],stdout=OUTPUT)
    OUTPUT.close()
    #call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/CTCF"+"/CTCF_In_Insulated_Neighbourhoods.tss",args.name.rstrip("/")+"/CTCF"+"/CTCF_In_Insulated_Neighbourhoods.png","TADs","0",str(args.ylimit)])
    #call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/CTCF"+"/CTCF_Outside_Insulated_Neighbourhoods.tss",args.name.rstrip("/")+"/CTCF"+"/CTCF_Outside_Insulated_Neighbourhoods.png",
    #     "NonTADs","0",str(args.ylimit)])
    #call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/CTCF"+"/CTCF_In_Insulated_Neighbourhoods.tss",args.name.rstrip("/")+"/CTCF"+"/CTCF_Outside_Insulated_Neighbourhoods.tss",
    #      args.name.rstrip("/")+"/CTCF"+"/CTCF_TADs.png","CTCF  sites in TAD boundaries","CTCF  sites outside TAD boundaries","0",str(args.ylimit)])
#########################################################################
def tf_gtrd_1000sites(args,avg_coverage):
    print("Analyze Transcription factors GTRD")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only_1000sites"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only_1000sites")
    target_list = glob.glob("./Ref/GTRD_1000sites/*.bed")
    for tf in target_list:
        tf_name = os.path.basename(tf[:-4])
        if os.path.isfile(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only_1000sites"+"/"+tf_name+".tss"):
            print("Skip "+tf_name)
            continue
        OUTPUT = open(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only_1000sites"+"/"+tf_name+".tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-m","100000","-limit","30","-bed",tf,"-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only_1000sites"+"/"+tf_name+".tss",args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only_1000sites"+"/"+tf_name+".png",tf_name,"0",str(args.ylimit)])
        #call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".tss","./Ref/TranscriptionFactors/MergedMaleProfiles/"+tf+".tss",
        #    args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP/"+tf+"_control.png",tf+" ("+args.name+")",tf+" (MergedMale)","0",str(args.ylimit)])
#########################################################################
def tf_gtrd(args,avg_coverage):
    print "Analyze Transcription factors GTRD"
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only")
    if args.hg38:
        target_list = glob.glob("./Ref/GTRD/hg38/*hg38.bed")
    else:
        target_list = glob.glob("./Ref/GTRD/*hg19.bed")
    for tf in target_list:
        tf_name = os.path.basename(tf[:-4])
        if os.path.isfile(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only"+"/"+tf_name+".tss"):
            print "Skip "+tf_name
            continue
        OUTPUT = open(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only"+"/"+tf_name+".tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-m","100000","-limit","30","-bed",tf,"-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only"+"/"+tf_name+".tss",args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only"+"/"+tf_name+".png",tf_name,"0",str(args.ylimit)])
        #call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".tss","./Ref/TranscriptionFactors/MergedMaleProfiles/"+tf+".tss",
        #    args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP/"+tf+"_control.png",tf+" ("+args.name+")",tf+" (MergedMale)","0",str(args.ylimit)])

#########################################################################
# TSS
def tss(args,avg_coverage):
    print "Analyze HK vs. Unexpr. TSS"
    if not os.path.isdir(args.name.rstrip("/")+"/TSS"):
        os.mkdir(args.name.rstrip("/")+"/TSS")
    OUTPUT = open(args.name.rstrip("/")+"/TSS"+"/HK_APPRIS_isoforms.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TSS/Housekeeping_APPRIS_isos_hg19_positions.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/TSS"+"/HK.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TSS/Housekeeping.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/TSS"+"/FANTOM_lower01.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TSS/Fantomlower01.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/TSS/HK.tss",args.name.rstrip("/")+"/TSS/FANTOM_lower01.tss",
            args.name.rstrip("/")+"/TSS/HK_vs_Unexpr.png","Housekeeping TSS","Unexpressed TSS","0",str(args.ylimit)])

####################################################################################################
# AndrogenReceptor
def androgen(args,avg_coverage):
    print "Analyze Androgen Receptor Binding sites"
    if not os.path.isdir(args.name.rstrip("/")+"/AR"):
        os.mkdir(args.name.rstrip("/")+"/AR")
    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_TARBS_All.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/TARBS.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_NARBS_All.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/NARBS.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_TARBS_GTRD_All.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/AR_GTRD_TARBS_intersect.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_NARBS_GTRD_All.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/AR_GTRD_NARBS_intersect.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_TARBS_GTRD_50perc.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/AR_GTRD_50perc_TARBS_intersect.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_NARBS_GTRD_50perc.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/AR_GTRD_50perc_NARBS_intersect.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/AR/AR_TARBS_All.tss",args.name.rstrip("/")+"/AR/AR_NARBS_All.tss",
            args.name.rstrip("/")+"/TSS/TARBS_vs_ARBS.png","T-AR binding sites","N-AR binding sites","0",str(args.ylimit)])
    call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/AR/AR_TARBS_GTRD_All.tss",args.name.rstrip("/")+"/AR/AR_TARBS_GTRD_All.tss",
            args.name.rstrip("/")+"/TSS/TARBS_vs_ARBS_GTRDintersect.png","T-AR binding sites (GTRD intersect)","N-AR binding sites (GTRD intersect)","0",str(args.ylimit)])
    call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/AR/AR_TARBS_GTRD_50perc.tss",args.name.rstrip("/")+"/AR/AR_NARBS_GTRD_50perc.tss",
            args.name.rstrip("/")+"/TSS/TARBS_vs_ARBS_GTRDintersect_50perc.png","T-AR binding sites (GTRD intersect,50perc)","N-AR binding sites (GTRD intersect,50perc)","0",str(args.ylimit)])


####################################################################################################
# Check for binding sites proximal and distal to Transcription start sites
def tf_tss(args,avg_coverage):
    print "Analyze distal and proximal TF binding sites"
    if not os.path.isdir(args.name.rstrip("/")+"/TSS_TF"):
        os.mkdir(args.name.rstrip("/")+"/TSS_TF")

    target_list = glob.glob("./Ref/TranscriptionFactors/TSS_intersects/*distal.bed")
    for tf in target_list:
        tf_name = os.path.basename(tf[:-10])
        proximal_tf = tf[:-10]+"proximal.bed"
        if os.path.isfile(args.name.rstrip("/")+"/TSS_TF/"+tf_name+"distal.tss"):
            print "Skip "+tf_name
            continue
        OUTPUT = open(args.name.rstrip("/")+"/TSS_TF/"+tf_name+"distal.tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed",tf,"-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
        OUTPUT.close()
        if os.path.isfile(args.name.rstrip("/")+"/TSS_TF/"+tf_name+"proximal.tss"):
            print "Skip "+tf_name
            continue
        OUTPUT = open(args.name.rstrip("/")+"/TSS_TF/"+tf_name+"proximal.tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed",proximal_tf,"-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/TSS_TF/"+tf_name+"distal.tss",args.name.rstrip("/")+"/TSS_TF/"+tf_name+"proximal.tss",
            args.name.rstrip("/")+"/TSS_TF/"+tf_name+".png",tf_name+" distal to TSS (>2kbp)",tf_name+" proximal to TSS (<2kbp)","0",str(args.ylimit)])

####################################################################################################

if args.analysis == "all":
    ctcf(args,avg_coverage)
    tf_gtrd_chip_only(args,avg_coverage)
    tss(args,avg_coverage)
elif args.analysis == "tss":
    tss(args,avg_coverage)
elif args.analysis == "androgen":
    androgen(args,avg_coverage)
elif args.analysis == "ctcf":
    ctcf(args,avg_coverage)
elif args.analysis == "tf_gtrd":
    tf_gtrd(args,avg_coverage)
elif args.analysis == "tf_gtrd_1000sites":
    tf_gtrd_1000sites(args,avg_coverage)
else:
    print "Unknown analysis type"
    print " Use any of:"
    print "   -) all"
    print "   -) ctcf"
    print "   -) androgen"
    print "   -) tf_gtrd"
    print "   -) tf_gtrd_1000sites"
    print "   -) tf_tss"

    
