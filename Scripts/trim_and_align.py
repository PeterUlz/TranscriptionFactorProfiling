#! /usr/bin/env python

# Prepare FastQ for nucleosome analyses
#  -) trim FastQ files
#  -) align trimmed reads
#  -) PCR dupliacte removal


#version 0.1: First draft
#version 0.2: Show step information in help output

import subprocess
import sys
import argparse
import os
import time
import shutil

version="0.2"
# Check command line arguments ###################################################################################
def checkArgs(args):
    if not os.path.isfile(args.fastq_file):
        errorExit("FastQ file not found")

# Error and Exit ###################################################################################
def errorExit(message):
   print "Error"
   print message
   print "Exit"
   sys.exit(1)

# Show steps ###################################################################################
def showSteps():
   epilog="-) Step1 create directory and create MD5 file of input\n"
   epilog += "-) Step2 Trim fastq\n"
   epilog += "-) Step3 alignment and conversion to BAM\n"
   epilog += "-) Step4 remove PCR duplicates\n"
   return epilog
   
# Logging module ###################################################################################
def logging(HANDLE, message):
    HANDLE.write(time.strftime("%d/%m/%Y:  %H:%M:%S  : "+message+"\n"))
    print(time.strftime("%d/%m/%Y:  %H:%M:%S  : "+message))
    HANDLE.flush()

# Error logging module ###################################################################################
def error_logging(HANDLE, message):
    HANDLE.write(time.strftime("%d/%m/%Y:  %H:%M:%S  : "+message+"\n"))
    print(time.strftime("%d/%m/%Y:  %H:%M:%S  : "+message))
    HANDLE.flush()
    sys.exit(1)


########################################################################################################
# Step1 create directory and create MD5 file of input
def step1(args,proj_dir):
    return_val=subprocess.call(["mkdir",proj_dir])
    if return_val != 0:
        sys.stderr.write("Cannot create directory "+proj_dir+"\n")
        sys.exit(1)



    MD5_OUTPUT = open(proj_dir +"/"+args.name+".md5","w")
    subprocess.call(["md5sum",args.fastq_file],stdout=MD5_OUTPUT)
    MD5_OUTPUT.close()


########################################################################################################
# Initialize analysis files and set up logging
def initalize_logging(args,proj_dir,version):
    ERR_LOG = open(proj_dir+"/"+args.name+".err.log","w")
    OUT_LOG = open(proj_dir+"/"+args.name+".out.log","w")
    logging(OUT_LOG,"Analysis Pipeline version: "+version)
    logging(OUT_LOG,"Starting Analysis")
    logging(OUT_LOG,"Fastq file: "+args.fastq_file)
    logging(OUT_LOG,"Sample name file: "+args.name)
    logging(OUT_LOG,"Output Directory: "+args.outdir)
    logging(OUT_LOG,"Temporary Directory: "+args.temp_dir)

    if args.keep:
        logging(OUT_LOG,"Keeping temporary files")
    #logging(OUT_LOG,"Step1 create directory and create MD5 file of input")

    return ERR_LOG,OUT_LOG
########################################################################################################
# Step 2 Trim fastq
def step2(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames):
    logging(OUT_LOG,"Step2 Read trimming")

    unzipped_fastq=subprocess.Popen(["zcat",args.fastq_file],stdout=subprocess.PIPE)
    unzip=subprocess.call([script_dir+"/Software/fastx_trimmer","-Q33","-f","53","-l","113","-z","-o",filenames["trimmed_fastq"]],stdin=unzipped_fastq.stdout)
    unzipped_fastq.wait()

########################################################################################################
# Step3 alignment and conversion to BAM
def step3(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames):
    logging(OUT_LOG,"Step3   a) alignment")
    SAM = open(filenames["sam_file"],"w")
    aln_value=subprocess.call([script_dir+"/Software/bwa","mem","-M","-t",str(args.threads),script_dir+"/Ref/hg19",filenames["trimmed_fastq"]],stderr=ERR_LOG,stdout=SAM)
    SAM.close()
    if aln_value != 0:
        error_logging(ERR_LOG,"Error in alignment")

    if not args.keep:
        subprocess.call(["rm",filenames["trimmed_fastq"]])


    logging(OUT_LOG,"Step3   b) conversion to BAM")
    bam_conversion=subprocess.call(["java","-jar",script_dir+"/Software/picard.jar","SortSam","SO=coordinate","INPUT="+filenames["sam_file"],"OUTPUT="+filenames["bam_file"],"VALIDATION_STRINGENCY=LENIENT","CREATE_INDEX=true"],stderr=ERR_LOG,stdout=OUT_LOG)
    if bam_conversion != 0:
        error_logging(ERR_LOG,"Error in BAM conversion")
    if not args.keep:
        subprocess.call(["rm",filenames["sam_file"]])

    logging(OUT_LOG,"Step3   c) BAM stats")
    STATS =open(filenames["bam_file_rmdup_stats"],"w")
    bam_stats=subprocess.call([script_dir+"/Software/samtools","stats",filenames["bam_file"]],stderr=ERR_LOG,stdout=STATS)
    if bam_stats != 0:
        error_logging(ERR_LOG,"Error in BAM stats")
    

######################################################################################################## 
# Step4 remove PCR duplicates
def step4(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames):
    logging(OUT_LOG,"Step4 Remove duplicates and index")
    rmdup_value=subprocess.call([script_dir+"/Software/samtools","rmdup","-s",filenames["bam_file"],filenames["bam_file_rmdup"]],stderr=ERR_LOG,stdout=OUT_LOG)
    if rmdup_value != 0:
        error_logging(ERR_LOG,"Error in PCR duplicate removal")

    if not args.keep:
        subprocess.call(["rm",filenames["bam_file"]])

    rmdup_value=subprocess.call([script_dir+"/Software/samtools","index",filenames["bam_file_rmdup"]],stderr=ERR_LOG,stdout=OUT_LOG)
    if rmdup_value != 0:
        error_logging(ERR_LOG,"Error in BAM indexing")



######################################################################################################## 
# Step5 analyze TSS profile in housekeeping vs. unexpressed genes and plot in R
def step5(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames):
    logging(OUT_LOG,"Step5 Calculate TSS profiles of housekeeping genes vs. unexpressed genes")
    logging(OUT_LOG,"Step5    a) Housekeeping genes")
    HOUSEKEEPING=open(filenames["housekeeping_profile"],"w")
    tss_profile=subprocess.call([script_dir+"/Scripts/analyze_TSS_coverage.py","-rg",script_dir+"/Ref/refSeq_extended_names_strand.bed","-norm",
          "-gl",script_dir+"/Ref/GeneLists/HK_gene_names.txt","-m","0","-b",filenames["bam_file_rmdup"],"-t",str(args.threads),"-tmp",args.temp_dir],stderr=ERR_LOG,stdout=HOUSEKEEPING)
    HOUSEKEEPING.close()
    if tss_profile != 0:
        error_logging(ERR_LOG,"Error in TSS profiling")


    logging(OUT_LOG,"Step5    b) Unexpressed genes")
    UNEXPRESSED=open(filenames["unexpressed_profile"],"w")
    tss_profile=subprocess.call([script_dir+"/Scripts/analyze_TSS_coverage.py","-rg",script_dir+"/Ref/refSeq_extended_names_strand.bed","-norm",
          "-gl",script_dir+"/Ref/GeneLists/Fantom5_all_lower0.1.txt","-m","0","-b",filenames["bam_file_rmdup"],"-t",str(args.threads),"-tmp",args.temp_dir],stderr=ERR_LOG,stdout=UNEXPRESSED)
    UNEXPRESSED.close()
    if tss_profile != 0:
        error_logging(ERR_LOG,"Error in TSS profiling")



    logging(OUT_LOG,"Step5    c) Plotting")
    INFILE1 = open(script_dir+"/Scripts/plot_TSS_profile.R","r")
    subprocess.call(["R","--slave","--args",filenames["housekeeping_profile"],filenames["unexpressed_profile"],filenames["profile_plot"]],stdin=INFILE1,stderr=ERR_LOG,stdout=OUT_LOG)
    INFILE1.close()


######################################################################################################## 
# Step6 extract coverage parameters for expression prediction
def step6(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames):
    logging(OUT_LOG,"Step6 Extract coverage parameters")
    logging(OUT_LOG,"Step6    a) 2K-TSS Coverage")
    TSS_2K =open(filenames["coverage_2k_tss"],"w")
    subprocess.call([script_dir+"/Scripts/analyze_all_genes_by_TSS_coverage_norm_LOG2.py","-rg",script_dir+"/Ref/refSeq_extended_names_strand.bed",
                   "-b",filenames["bam_file_rmdup"],"-t",str(args.threads),"-gl",script_dir+"/Ref/GeneLists/AllGenes_NMonly.txt","-norm","-norm-file",args.cna_file,"-tmp","/tmp/"],stdout=TSS_2K,stderr=ERR_LOG)
    TSS_2K.close()

    logging(OUT_LOG,"Step6    b) NDR Coverage")
    TSS_NDR =open(filenames["coverage_ndr_tss"],"w")
    subprocess.call([script_dir+"/Scripts/analyze_all_genes_by_TSS_coverage_norm_LOG2.py","-s","150","-e","50","-rg",script_dir+"/Ref/refSeq_extended_names_strand.bed",
               "-b",filenames["bam_file_rmdup"],"-t",str(args.threads),"-gl",script_dir+"/Ref/GeneLists/AllGenes_NMonly.txt","-norm","-norm-file",args.cna_file,"-tmp","/tmp/"],stdout=TSS_NDR,stderr=ERR_LOG)
    TSS_NDR.close()

######################################################################################################## 
# Step7 expression prediction
def step7(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames):
    logging(OUT_LOG,"Step7 Expression prediction")
    return_val=subprocess.call(["mkdir",filenames["prediction_folder"]])   
    INFILE_PRED = open(script_dir+"/Scripts/prediction_svm_housekeeping.R","r") 
    r_return_value=subprocess.call(["R","--slave","--args",filenames["coverage_2k_tss"],filenames["coverage_ndr_tss"],filenames["prediction_folder"],script_dir+"/Ref"],stdin=INFILE_PRED,stderr=ERR_LOG,stdout=OUT_LOG)
    if r_return_value != 0:
        error_logging(ERR_LOG,"Error in expression prediction")

    INFILE_PRED.close()

   

###################################################################################
script_dir = os.path.dirname(os.path.realpath(__file__))
# Parse command line arguments ###################################################################################
epilog_text=showSteps()
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription end',epilog=epilog_text,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-fq','--fastq', dest='fastq_file', 
                   help='Fastq file (for NextSeq: Specify Lane 1 Fastq File)',required=True)
parser.add_argument('-s','--sample-name', dest='name',
                   help='Sample name to be used subsequently',required=True)
parser.add_argument('-g','--gender', dest='gender',
                   help='Gender of the sample [either m or f]',required=True, choices=["m","f"])
parser.add_argument('-o','--out-dir', dest='outdir',
                   help='Output Directory [default .]',default=".")
parser.add_argument('-k','--keep-temp', dest='keep',
                   help='Keep temporary files',action="store_true")
parser.add_argument('-t','--threads', dest='threads',
                   help='No. threads for alignment [default: 1]',type=int,default=1)
parser.add_argument('-cna','--cna-file', dest='cna_file',
                   help='File containing log2-ratios of copy-number segments (.segments file)',required=True)
parser.add_argument('-tmp','--temp-dir', dest='temp_dir',
                   help='Temporary Directory',default="/tmp/")
parser.add_argument('-step','--start-step', dest='start_step',
                   help='Start at this step',type=int,default=1)

args = parser.parse_args()
checkArgs(args)
if args.outdir[-1] != "/":
    proj_dir = args.outdir+"/"+args.name
else:
    proj_dir = args.outdir+args.name

#####################################################################################################
# Define filenames
filenames = dict()
filenames["trimmed_fastq"]=proj_dir+"/"+args.name+".trimmed.fastq.gz"
filenames["sam_file"] = proj_dir+"/"+args.name+".sam"
filenames["bam_file"] = proj_dir+"/"+args.name+".bam"
filenames["bam_file_rmdup"] = proj_dir+"/"+args.name+".rmdup.bam"
filenames["bam_file_rmdup_stats"] = proj_dir+"/"+args.name+".rmdup.bam.stats"

filenames["housekeeping_profile"] = proj_dir+"/"+args.name+".housekeeping.tss"
filenames["unexpressed_profile"] = proj_dir+"/"+args.name+".unexpressed.tss"
filenames["profile_plot"] = proj_dir+"/"+args.name+".tss_profile.png"

filenames["coverage_2k_tss"] = proj_dir+"/"+args.name+".2k-tss.txt"
filenames["coverage_ndr_tss"] = proj_dir+"/"+args.name+".ndr-tss.txt"

filenames["prediction_folder"] = proj_dir+"/Prediction"


#####################################################################################################
# Start analysis
if (args.start_step == 1):
    step1(args,proj_dir)

ERR_LOG,OUT_LOG=initalize_logging(args,proj_dir,version)
if (args.start_step > 1):
    logging(OUT_LOG, "Skip creating directory and MD5 sum")

if (args.start_step <= 2):
    step2(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames)
else:
    logging(OUT_LOG, "Skip read trimming")

if (args.start_step <= 3):
    step3(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames)
else:
    logging(OUT_LOG, "Skip alignment and BAM conversion")

if (args.start_step <= 4):
    step4(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames)
else:
    logging(OUT_LOG, "Skip PCR duplicate removal")

if (args.start_step <= 5):
    step5(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames)
else:
    logging(OUT_LOG, "Skip TSS profile generation")

if (args.start_step <= 6):
    step6(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames)
else:
    logging(OUT_LOG, "Skip parameter extraction")

if (args.start_step <= 7):
    step7(ERR_LOG,OUT_LOG,args,proj_dir,script_dir,filenames)
else:
    logging(OUT_LOG, "Skip expression prediction")

