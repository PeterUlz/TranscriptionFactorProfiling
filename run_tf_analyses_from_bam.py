#! /usr/bin/python

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
# TODO:
#  -) better implement hg38 support
#  -) add standard TSS (housekeeping and unexpressed) analysis
####################################################################################################
# get genomewide coverage from bedtools genomecoverage
if args.calc_cov:
    print "Calc avg. coverage"
    OUTPUT=open(args.name.rstrip("/")+"/"+args.name+".coverage","w")
    if args.hg38:
        call(["bedtools","genomecov","-ibam",args.bam_file,"-g","./Ref/hg38.chrom_sizes.txt"],stdout=OUTPUT)
    else:
        call(["bedtools","genomecov","-ibam",args.bam_file,"-g","./Ref/hg19.chrom_sizes.txt"],stdout=OUTPUT)
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

####################################################################################################
# G Quadruplexes
def quadruplex(args,avg_coverage):
    print "Analyze G-quadruplexes"
    if not os.path.isdir(args.name.rstrip("/")+"/G-quadruplexes"):
        os.mkdir(args.name.rstrip("/")+"/G-quadruplexes")
    OUTPUT = open(args.name.rstrip("/")+"/G-quadruplexes"+"/4Tetrads.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/G-quadruplexes/Quadruplex_4_tetrads_both_strands.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/G-quadruplexes"+"/3Tetrads.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/G-quadruplexes/Quadruplex_3_tetrads_both_strands.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-m","100000","-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/G-quadruplexes"+"/5Tetrads.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/G-quadruplexes/Quadruplex_5_tetrads_both_strands.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage_3sample.R",args.name.rstrip("/")+"/G-quadruplexes"+"/3Tetrads.tss",args.name.rstrip("/")+"/G-quadruplexes"+"/4Tetrads.tss",args.name.rstrip("/")+"/G-quadruplexes"+"/5Tetrads.tss",args.name.rstrip("/")+"/G-quadruplexes"+"/G-quadruplexes.png","0",str(args.ylimit)])


####################################################################################################
# PolyAsites
def polya(args,avg_coverage):
    print "Analyze PolyA sites"
    if not os.path.isdir(args.name.rstrip("/")+"/PolyA"):
        os.mkdir(args.name.rstrip("/")+"/PolyA")
    OUTPUT = open(args.name.rstrip("/")+"/PolyA"+"/apaDb_all.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/PolyA-Sites/apaDb.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/PolyA"+"/apaDb_all.tss",args.name.rstrip("/")+"/PolyA"+"/apaDb_all.png","All PolyA sites","0",str(args.ylimit)])
    OUTPUT = open(args.name.rstrip("/")+"/PolyA"+"/apaDb_over100reads.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/PolyA-Sites/apaDb_over100reads.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/PolyA"+"/apaDb_over100reads.tss",args.name.rstrip("/")+"/PolyA"+"/apaDb_over100reads.png","All PolyA sites >100 MACE reads","0",str(args.ylimit)])
    OUTPUT = open(args.name.rstrip("/")+"/PolyA"+"/apaDb_HK.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/PolyA-Sites/apaDb_HK.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/PolyA"+"/apaDb_HK.tss",args.name.rstrip("/")+"/PolyA"+"/apaDb_HK.png","All PolyA sites in Housekeeping genes","0",str(args.ylimit)])
    OUTPUT = open(args.name.rstrip("/")+"/PolyA"+"/apaDb_HK_over100reads.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/PolyA-Sites/apaDb_HK_over100reads.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/PolyA"+"/apaDb_HK.tss",args.name.rstrip("/")+"/PolyA"+"/apaDb_HK_over100reads.png","All PolyA sites in Housekeepng genes >100 MACE reads","0",str(args.ylimit)])
    call(["Rscript","./Scripts/plot_MotifCoverage_3sample.R",args.name.rstrip("/")+"/PolyA"+"/apaDb_all.tss",args.name.rstrip("/")+"/PolyA"+"/apaDb_over100reads.tss",args.name.rstrip("/")+"/PolyA"+"/apaDb_HK_over100reads.tss",args.name.rstrip("/")+"/PolyA"+"/PolyA-sites.png","0",str(args.ylimit)])
###################################################################################################
def tf_gtrd_1000sites(args,avg_coverage):
    print("Analyze Transcription factors GTRD")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only_1000sites"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only_1000sites")
    target_list = glob.glob("./Ref/TranscriptionFactors/GTRD_1000sites/*.bed")
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

####################################################################################################
# X inactivation profiles around transcription start sites
def x_inactivation(args,avg_coverage):
    print "Analyze X_inactivation"
    if not os.path.isdir(args.name.rstrip("/")+"/X_inactivation"):
        os.mkdir(args.name.rstrip("/")+"/X_inactivation")
    #x_coverage = getChromCoverage("chrX",args)
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_all.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_allGenes.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_all.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_all.png","All TSS on chrX","0",str(args.ylimit)])
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_HK.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_Housekeeping_TSS.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_HK.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_HK.png","HK TSS on chrX","0",str(args.ylimit)])
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_FPKM_over8.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_FPKM_over8.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_FPKM_over8.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_FPKM_over8.png","TSS in Genes FPKM >8 chrX","0",str(args.ylimit)])
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_FPKM_lower01.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_FPKM_lower01.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_FPKM_lower01.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_FPKM_lower01.png","TSS in Genes FPKM >8 chrX","0",str(args.ylimit)])
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_Escapee_genes.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_Escapee_genes.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_Escapee_genes.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_Escapee_genes.png","TSS in genes that escape XCI","0",str(args.ylimit)])
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_Escapee_and_Mostly_genes.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_Escapee_and_Mostly_genes.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_Escapee_and_Mostly_genes.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_Escapee_and_Mostly_genes.png","TSS in genes that mostly escape XCI","0",str(args.ylimit)])
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_sensitive.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_XCI_sensitive.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_sensitive.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_sensitive.png","TSS in genes that are sensitive to XCI","0",str(args.ylimit)])
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_sensitive_and_mostly.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_XCI_sensitive_and_mostly.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_sensitive_and_mostly.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_sensitive_and_mostly.png","TSS in genes that are mostly sensitive to XCI","0",str(args.ylimit)])
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_mostly_escape_FPKMover1.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_XCI_mostly_escape_FPKMover1.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_mostly_escape_FPKMover1.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_mostly_escape_FPKMover1.png","TSS in genes that mostly escape XCI with FPKM > 1","0",str(args.ylimit)])
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_mostly_sensitive_FPKMover1.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_XCI_mostly_sensitive_FPKMover1.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_mostly_sensitive_FPKMover1.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_mostly_sensitive_FPKMover1.png","TSS in genes that are mostly sensitive to XCI with FPKM > 1","0",str(args.ylimit)])
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_escape_Tukiainen.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_XCI_escape_Tukiainen.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_escape_Tukiainen.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_escape_Tukiainen.png","TSS in genes that escape XCI","0",str(args.ylimit)])
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_inactive_Tukiainen.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_XCI_inactive_Tukiainen.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_inactive_Tukiainen.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_inactive_Tukiainen.png","TSS of inactive genes","0",str(args.ylimit)])
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_variable_Tukiainen.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-limit","50","-bed","./Ref/X_inactivation/chromX_XCI_variable_Tukiainen.bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_variable_Tukiainen.tss",args.name.rstrip("/")+"/X_inactivation"+"/chromX_XCI_variable_Tukiainen.png","TSS of variable genes","0",str(args.ylimit)])

####################################################################################################
# Transcription factors Clasical
def tf(args,avg_coverage):
    print "Analyze Transcription factors"
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP")
    tf_list = ['MYB','ATF2','CEBPB','MYC','TP53','E2F4','MAFK','BRCA1','SMAD1','NFYA','ETV6','AR','JUND','MAX','NFAC1','PBX3','RELB','MXI1','RFX5','SPI1','TCF7','SRF','TBX21','STAT1']
    for tf in tf_list:
        OUTPUT = open(args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/ENCODE_ChIP/"+tf+".txt","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".tss",args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".png",tf,"0",str(args.ylimit)])
        call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".tss","./Ref/TranscriptionFactors/MergedMaleProfiles/"+tf+".tss",
            args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP/"+tf+"_control.png",tf+" ("+args.name+")",tf+" (MergedMale)","0",str(args.ylimit)])

####################################################################################################
# Transcription factors Clasical (read Paper by Nuria Lopez-Bigas)
def exonIntron(args,avg_coverage):
    print "Analyze Exon Intron Boundaries"
    if not os.path.isdir(args.name.rstrip("/")+"/ExonIntronBoundaries"):
        os.mkdir(args.name.rstrip("/")+"/ExonIntronBoundaries")
    OUTPUT = open(args.name.rstrip("/")+"/ExonIntronBoundaries"+"/All_exons.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/ExonIntronBoundaries/exons_noXY.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/ExonIntronBoundaries"+"/All_exons.tss",args.name.rstrip("/")+"/ExonIntronBoundaries"+"/All_exons.png","Exon-Intron boundary","0",str(args.ylimit)])

###################################################################################################
def tf_gtrd_fimo(args,avg_coverage):
    print "Analyze Transcription factors GTRD"
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD")
    target_list = glob.glob("./Ref/TranscriptionFactors/GTRD_FIMO_intersect/*.txt")
    for tf in target_list:
        tf_name = os.path.basename(tf[:-4])
        if os.path.isfile(args.name.rstrip("/")+"/TranscriptionFactors/GTRD"+"/"+tf_name+".tss"):
            print "Skip "+tf_name
            continue
        OUTPUT = open(args.name.rstrip("/")+"/TranscriptionFactors/GTRD"+"/"+tf_name+".tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed",tf,"-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/TranscriptionFactors/GTRD"+"/"+tf_name+".tss",args.name.rstrip("/")+"/TranscriptionFactors/GTRD"+"/"+tf_name+".png",tf_name,"0",str(args.ylimit)])
        #call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".tss","./Ref/TranscriptionFactors/MergedMaleProfiles/"+tf+".tss",
        #    args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP/"+tf+"_control.png",tf+" ("+args.name+")",tf+" (MergedMale)","0",str(args.ylimit)])
###################################################################################################
def tf_gtrd_chip_only(args,avg_coverage):
    print "Analyze Transcription factors GTRD"
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_ChIP_Only")
    if args.hg38:
        target_list = glob.glob("./Ref/TranscriptionFactors/GTRD/hg38/*hg38.bed")
    else:
        target_list = glob.glob("./Ref/TranscriptionFactors/GTRD/*hg19.bed")
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
###################################################################################################
def tf_tissue(args,avg_coverage):
    print "Analyze Transcription factors from Virtual ChIP-seq"
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors/VirChIP"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors/VirChIP")
    tf_select = args.tf
    print "  -) search for "+tf_select+" ref files"
    target_list = glob.glob("./Ref/TranscriptionFactors/bedPredictions/hg19/hg19_GTRD_intersect/*"+tf_select+"*hg19.bed.GTRDintersect.bed")
    for tf in target_list:
        tf_name = os.path.basename(tf[:-8])
        if os.path.isfile(args.name.rstrip("/")+"/TranscriptionFactors/VirChIP/"+tf_name+".GTRD_interect.tss"):
            print "Skip "+tf_name
            continue
        OUTPUT = open(args.name.rstrip("/")+"/TranscriptionFactors/VirChIP"+"/"+tf_name+".tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-m","100000","-limit","30","-bed",tf,"-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/TranscriptionFactors/VirChIP/"+tf_name+".tss",args.name.rstrip("/")+"/TranscriptionFactors/VirChIP"+"/"+tf_name+".png",tf_name,"0",str(args.ylimit)])
        #call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".tss","./Ref/TranscriptionFactors/MergedMaleProfiles/"+tf+".tss",
        #    args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP/"+tf+"_control.png",tf+" ("+args.name+")",tf+" (MergedMale)","0",str(args.ylimit)])

###################################################################################################
def lola(args,avg_coverage):
    print "Analyze epigenetic regions from LOLA database"
    if not os.path.isdir(args.name.rstrip("/")+"/LOLA"):
        os.mkdir(args.name.rstrip("/")+"/LOLA")
    target_list = glob.glob("./Ref/LOLA/LOLACore/hg19/*/regions/*.bed")
    for tf in target_list:
        tf_name = os.path.basename(tf[:-4])
        if os.path.isfile(args.name.rstrip("/")+"/LOLA/"+tf_name+".tss"):
            print "Skip "+tf_name
            continue
        OUTPUT = open(args.name.rstrip("/")+"/LOLA/"+tf_name+".tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-m","100000","-limit","30","-bed",tf,"-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/LOLA/"+tf_name+".tss",args.name.rstrip("/")+"/LOLA/"+tf_name+".png",tf_name,"0",str(args.ylimit)])
        #call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".tss","./Ref/TranscriptionFactors/MergedMaleProfiles/"+tf+".tss",
        #    args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP/"+tf+"_control.png",tf+" ("+args.name+")",tf+" (MergedMale)","0",str(args.ylimit)])


###################################################################################################
def tf_gtrd_fimo_q30(args,avg_coverage):
    print "Analyze Transcription factors GTRD"
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Q30"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Q30")
    target_list = glob.glob("./Ref/TranscriptionFactors/GTRD_FIMO_intersect/*.txt")
    for tf in target_list:
        tf_name = os.path.basename(tf[:-4])
        if os.path.isfile(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Q30"+"/"+tf_name+".tss"):
            print "Skip "+tf_name
            continue
        OUTPUT = open(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Q30"+"/"+tf_name+".tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed",tf,"-cov",str(avg_coverage),"-norm","-mapq","30","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Q30"+"/"+tf_name+".tss",args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Q30"+"/"+tf_name+".png",tf_name,"0",str(args.ylimit)])
        #call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".tss","./Ref/TranscriptionFactors/MergedMaleProfiles/"+tf+".tss",
        #    args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP/"+tf+"_control.png",tf+" ("+args.name+")",tf+" (MergedMale)","0",str(args.ylimit)])
###################################################################################################
def tf_gtrd_fimo_heatmap(args,avg_coverage):
    print "Analyze Transcription factors GTRD"
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors")
    if not os.path.isdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Heatmap"):
        os.mkdir(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Heatmap")
    target_list = glob.glob("./Ref/TranscriptionFactors/GTRD_FIMO_intersect/*.txt")
    for tf in target_list:
        tf_name = os.path.basename(tf[:-4])
        if os.path.isfile(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Heatmap"+"/"+tf_name+".tss"):
            print "Skip "+tf_name
            continue
        OUTPUT = open(args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Heatmap"+"/"+tf_name+".tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam_heatmap.py","-m","100000","-bed",tf,"-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-n",args.name.rstrip("/")+"/TranscriptionFactors/GTRD_Heatmap"+"/"+tf_name+".png"],stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/TranscriptionFactors/GTRD"+"/"+tf_name+".tss",args.name.rstrip("/")+"/TranscriptionFactors/GTRD"+"/"+tf_name+".png",tf_name,"0",str(args.ylimit)])
        #call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/TranscriptionFactors/ENCODE_ChIP"+"/"+tf+".tss","./Ref/TranscriptionFactors/MergedMaleProfiles/"+tf+".tss",
 
####################################################################################################
# Histone Modifications:
def histone(args,avg_coverage):
    print "Analyze Histone Modifications"
    if not os.path.isdir(args.name.rstrip("/")+"/HistoneModifications"):
        os.mkdir(args.name.rstrip("/")+"/HistoneModifications")
    mod_list = ['H3K9ac','H3K79me2','H3K27ac','H3K4me2','H4K20me1','H3K36me3','H3K27me3','H3K4me3','H2AFZ','H3K4me1']
    for mod in mod_list:
        OUTPUT = open(args.name.rstrip("/")+"/HistoneModifications/"+mod+".tss","w")
        call(["./Scripts/analyze_coverage_around_position_log2_pysam_heatmap.py","-bed","./Ref/HistoneModifications/"+mod+".bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000","-n",args.name.rstrip("/")+"/HistoneModifications/"+mod+".heatmap.png"])
        call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/HistoneModifications/"+mod+".bed","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_MotifCoverage.R",args.name.rstrip("/")+"/HistoneModifications/"+mod+".tss",args.name.rstrip("/")+"/HistoneModifications/"+mod+".png",mod,"0",str(args.ylimit)])
####################################################################################################
# Enhancers:
def enhancer(args,avg_coverage):
    print "Analyze Enhancers"
    if not os.path.isdir(args.name.rstrip("/")+"/Enhancers"):
        os.mkdir(args.name.rstrip("/")+"/Enhancers")
    enhancer_list = ['basophil','b_cell','blood','brain','dendritic','hepatocyte','macrophage','monocyte','neutrophil','nk_cell','t_cell']
    
    for mod in enhancer_list:
        OUTPUT = open(args.name.rstrip("/")+"/Enhancers"+"/"+mod+".bed","w")
        call(["./Software/mosdepth","-t","4","-b","./Ref/Enhancers/slidebase_enhancers_"+mod+"_30perc_Expr_sort_noXY.bed",args.bam_file], stdout=OUTPUT)
        OUTPUT.close()
        OUTPUT = open(args.name.rstrip("/")+"/Enhancers"+"/"+mod+"_ref.bed","w")
        call(["./Software/mosdepth","-t","4","-b","./Ref/Enhancers/slidebase_enhancers_"+mod+"_30perc_Expr_sort_noXY_ref.bed",args.bam_file], stdout=OUTPUT)
        OUTPUT.close()
        call(["Rscript","./Scripts/plot_RegionCoverage_2sample.R",args.name.rstrip("/")+"/Enhancers"+"/"+mod+".bed",
              args.name.rstrip("/")+"/Enhancers"+"/"+mod+"_ref.bed",args.name.rstrip("/")+"/Enhancers"+"/"+mod,str(args.ylimit+1),str(avg_coverage)])

####################################################################################################
# Replication Origins:
def replication_origins(args,avg_coverage):
    print "Analyze Replication Origins"
    if not os.path.isdir(args.name.rstrip("/")+"/ReplicationOrigins"):
        os.mkdir(args.name.rstrip("/")+"/ReplicationOrigins")
    OUTPUT = open(args.name.rstrip("/")+"/ReplicationOrigins/RepliSeqPeaks.bed","w")
    call(["./Software/mosdepth","-t","4","-b","./Ref/RepliSeq/RepliSeqPeaks.bed",args.bam_file], stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/ReplicationOrigins/RepliSeqPeaks_ref.bed","w")
    call(["./Software/mosdepth","-t","4","-b","./Ref/RepliSeq/RepliSeqPeaks_random_ref.bed",args.bam_file], stdout=OUTPUT)
    OUTPUT.close()
    call(["Rscript","./Scripts/plot_RegionCoverage_2sample.R",args.name.rstrip("/")+"/ReplicationOrigins/RepliSeqPeaks.bed",
              args.name.rstrip("/")+"/ReplicationOrigins/RepliSeqPeaks_ref.bed",args.name.rstrip("/")+"/ReplicationOrigins/RepliSeq.png",str(args.ylimit+1),str(avg_coverage)])
####################################################################################################
# TSS gene sets:
def tss_gene_sets(args,avg_coverage):
    print "Analyze TSSs of gene sets"
    if not os.path.isdir(args.name.rstrip("/")+"/TSS_gene_sets"):
        os.mkdir(args.name.rstrip("/")+"/TSS_gene_sets")
    if not os.path.isdir(args.name.rstrip("/")+"/TSS_gene_sets/All"):
        os.mkdir(args.name.rstrip("/")+"/TSS_gene_sets/All")
    if not os.path.isdir(args.name.rstrip("/")+"/TSS_gene_sets/APPRIS"):
        os.mkdir(args.name.rstrip("/")+"/TSS_gene_sets/APPRIS")
    call(["./Scripts/TSS_gene_sets/analyze_gene_sets.py",args.bam_file,args.norm_log2,"./Ref/TSS_gene_sets/refSeq_extended_names_strand.bed",args.name.rstrip("/")+"/TSS_gene_sets/All/",str(avg_coverage)])
    call(["./Scripts/TSS_gene_sets/analyze_gene_sets.py",args.bam_file,args.norm_log2,"./Ref/TSS_gene_sets/ref_seq_genes_appris_hg19.bed",args.name.rstrip("/")+"/TSS_gene_sets/APPRIS/",str(avg_coverage)])


####################################################################################################
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
#    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_TARBS_All.tss","w")
#    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/TARBS.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
#    OUTPUT.close()
#    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_NARBS_All.tss","w")
#    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/NARBS.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
#    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_TARBS_GTRD_All.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/AR_GTRD_TARBS_intersect.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_NARBS_GTRD_All.tss","w")
    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/AR_GTRD_NARBS_intersect.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
    OUTPUT.close()
#    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_TARBS_GTRD_50perc.tss","w")
#    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/AR_GTRD_50perc_TARBS_intersect.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
#    OUTPUT.close()
#    OUTPUT = open(args.name.rstrip("/")+"/AR"+"/AR_NARBS_GTRD_50perc.tss","w")
#    call(["./Scripts/analyze_coverage_around_position_log2_pysam.py","-bed","./Ref/TranscriptionFactors/AndrogenReceptor/AR_GTRD_50perc_NARBS_intersect.bed","-m","100000","-cov",str(avg_coverage),"-norm","-norm-file",args.norm_log2,"-b",args.bam_file,"-s","1000","-e","1000"],stdout=OUTPUT)
#    OUTPUT.close()
    call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/AR/AR_TARBS_All.tss",args.name.rstrip("/")+"/AR/AR_NARBS_All.tss",
            args.name.rstrip("/")+"/TSS/TARBS_vs_ARBS.png","T-AR binding sites","N-AR binding sites","0",str(args.ylimit)])
    call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/AR/AR_TARBS_GTRD_All.tss",args.name.rstrip("/")+"/AR/AR_TARBS_GTRD_All.tss",
            args.name.rstrip("/")+"/TSS/TARBS_vs_ARBS_GTRDintersect.png","T-AR binding sites (GTRD intersect)","N-AR binding sites (GTRD intersect)","0",str(args.ylimit)])
    call(["Rscript","./Scripts/plot_MotifCoverage_2sample.R",args.name.rstrip("/")+"/AR/AR_TARBS_GTRD_50perc.tss",args.name.rstrip("/")+"/AR/AR_NARBS_GTRD_50perc.tss",
            args.name.rstrip("/")+"/TSS/TARBS_vs_ARBS_GTRDintersect_50perc.png","T-AR binding sites (GTRD intersect,50perc)","N-AR binding sites (GTRD intersect,50perc)","0",str(args.ylimit)])


####################################################################################################
# AndrogenReceptor
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
    quadruplex(args,avg_coverage)
    tf(args,avg_coverage)
    tf_gtrd_fimo(args,avg_coverage)
    tf_gtrd_fimo_q30(args,avg_coverage)
    tf_gtrd_chip_only(args,avg_coverage)
    polya(args,avg_coverage)
    x_inactivation(args,avg_coverage)
    histone(args,avg_coverage)
    enhancer(args,avg_coverage)
    replication_origins(args,avg_coverage)
    tss(args,avg_coverage)
elif args.analysis == "tss":
    tss(args,avg_coverage)
elif args.analysis == "androgen":
    androgen(args,avg_coverage)
elif args.analysis == "ctcf":
    ctcf(args,avg_coverage)
elif args.analysis == "quadruplex":
    quadruplex(args,avg_coverage)
elif args.analysis == "tf":
    tf(args,avg_coverage)
elif args.analysis == "tf_tss":
    tf_tss(args,avg_coverage)
elif args.analysis == "tf_gtrd":
    tf_gtrd_fimo(args,avg_coverage)
elif args.analysis == "tf_gtrd_q30":
    tf_gtrd_fimo_q30(args,avg_coverage)
elif args.analysis == "tf_gtrd_heatmap":
    tf_gtrd_fimo_heatmap(args,avg_coverage)
elif args.analysis == "tf_gtrd_chip_only":
    tf_gtrd_chip_only(args,avg_coverage)
elif args.analysis == "tf_gtrd_1000sites":
    tf_gtrd_1000sites(args,avg_coverage)
elif args.analysis == "tf_tissue":
    tf_tissue(args,avg_coverage)
elif args.analysis == "lola":
    lola(args,avg_coverage)
elif args.analysis == "polya":
    polya(args,avg_coverage)
elif args.analysis == "x_inactivation":
    x_inactivation(args,avg_coverage)
elif args.analysis == "exon_intron":
    exonIntron(args,avg_coverage)
elif args.analysis == "histone":
    histone(args,avg_coverage)
elif args.analysis == "enhancer":
    enhancer(args,avg_coverage)
elif args.analysis == "replication_origins":
    replication_origins(args,avg_coverage)
elif args.analysis == "tss_gene_sets":
    tss_gene_sets(args,avg_coverage)
else:
    print "Unknown analysis type"
    print " Use any of:"
    print "   -) all"
    print "   -) ctcf"
    print "   -) quadruplex"
    print "   -) tf"
    print "   -) histone"
    print "   -) enhancer"

    
