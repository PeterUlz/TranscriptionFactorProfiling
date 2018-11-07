#! /usr/bin/python

# Trim from paired-end BAM file
# this needs to be done in order to arrive at the "central 60bp" 
import pysam
import sys
import argparse

# Parse command line arguments ###################################################################################
parser = argparse.ArgumentParser(description='Analyze read depth in comparison to transcription start')
parser.add_argument('-i','--input-file', dest='in_file', 
                   help='Input BAM file',required=True)
parser.add_argument('-o','--output-file', dest='out_file',
                   help='Output BAM file',required=True)
parser.add_argument('-l','--length', dest='length',
                   help='Output length',type=int,required=True)
#parser.add_argument('-r','--ref-fasta', dest='ref',
#                   help='reference fasta file',required=True)
args = parser.parse_args()

bamF = pysam.AlignmentFile(args.in_file, 'rb')
bamOut = pysam.AlignmentFile(args.out_file, 'wb', template=bamF)
counter = 0

###########################################################################
for read in bamF.fetch():
    read_start = read.reference_start
    read_end = read.reference_end
    reverse = False
    if read.is_reverse:
        reverse = True 

    counter += 1
    if (counter  % 10000 == 0):
        sys.stdout.write("\r"+str(counter)+" reads analyzed")
        sys.stdout.flush()
    
    if not reverse:
        start = read_start + ((166-args.length) // 2)
    else:
        start = read_end - args.length - ((166-args.length) // 2)
        if start < 1:
            start = 1
    trimmed_read = pysam.AlignedSegment()
    trimmed_read.query_name = read.query_name
    trimmed_read.query_sequence="A"*args.length
    trimmed_read.flag = read.flag
    trimmed_read.reference_id = read.reference_id
    trimmed_read.reference_start = start
    trimmed_read.mapping_quality = read.mapping_quality
    trimmed_read.cigar = [(0,args.length)]
    #trimmed_read.next_reference_id = read.next_reference_id
    #trimmed_read.next_reference_start=read.next_reference_start
    trimmed_read.template_length=read.template_length
    trimmed_read.query_qualities = [37] * args.length
    trimmed_read.tags = read.tags
 
    bamOut.write(trimmed_read)

bamF.close()
bamOut.close()
pysam.sort("-o",args.out_file+".sorted.bam",args.out_file)
pysam.index(args.out_file+".sorted.bam")

