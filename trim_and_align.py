#! /usr/bin/python

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
    aln_value=subprocess.call([script_dir+"/Software/bwa","mem","-M","-t",str(args.threads),args.bwa_idx,filenames["trimmed_fastq"]],stderr=ERR_LOG,stdout=SAM)
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

###################################################################################
script_dir = os.path.dirname(os.path.realpath(__file__))
# Parse command line arguments ###################################################################################
epilog_text=showSteps()
parser = argparse.ArgumentParser(description='Trim and align cfDNA data',epilog=epilog_text,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-fq','--fastq', dest='fastq_file', 
                   help='Fastq file (for NextSeq: Specify Lane 1 Fastq File)',required=True)
parser.add_argument('-s','--sample-name', dest='name',
                   help='Sample name to be used subsequently',required=True)
parser.add_argument('-o','--out-dir', dest='outdir',
                   help='Output Directory [default .]',default=".")
parser.add_argument('-k','--keep-temp', dest='keep',
                   help='Keep temporary files',action="store_true")
parser.add_argument('-t','--threads', dest='threads',
                   help='No. threads for alignment [default: 1]',type=int,default=1)
parser.add_argument('-step','--start-step', dest='start_step',
                   help='Start at this step',type=int,default=1)
parser.add_argument('-idx','--bwa-index', dest='bwa_idx',
                   help='BWA index prefix',required=True)
parser.add_argument('-tmp','--tmp-dir', dest='temp_dir',
                   help='Directory for temporary file [default: /tmp]',default="/tmp")

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

