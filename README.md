# TranscriptionFactorProfiling
Profiling of transcription factor binding sites in cell-free DNA

Raw sequencing data is trimmed to contain bases 53-113 (i.e. the central 60bp in a hypothetical 166bp fragment), PCR duplicates are removed and trimmed read is aligned.

In a second step, coverage data around transcription factor binding sites are averaged and the activity for each transcription factor is calculated.

Transcription factor binding sites are taken from GTRD database

## Preparation:
Software needed:
* bedtools (v2.26.0)
* bwa (0.7.4-r385)
* samtools (v 1.3)
* fastx_trimmer
* picard 

Binaries should be put in ./Software/

## Trim and align
`
usage: trim_and_align.py [-h] -fq FASTQ_FILE -s NAME [-o OUTDIR] [-k]
                         [-t THREADS] [-step START_STEP]
`

## Calculate copy-number alterations
This is used in order to correct for local coverage differences due to copy-number alterations. 
The output of this is directly fed into the subsequent program, where coverage 
One option to do this is Plasma-Seq (https://github.com/PeterUlz/PlasmaSeq).
Other options could include: ichorCNA (https://github.com/broadinstitute/ichorCNA) or qDNAseq (https://bioconductor.org/packages/release/bioc/html/QDNAseq.html)

However, make sure the input file for the following steps looks like this:
<chr>TAB<start>TAB<end>TAB<Log2-ratio>
This corresponds to the .segments data produced by Plasma-Seq


## Run coverage analysis over all transcription factors in GTRD
`
usage: run_tf_analyses_from_bam.py [-h] -b BAM_FILE -o NAME  
                                   [-cov MEAN_COVERAGE] [-ylimit YLIMIT]  
                                   [-norm-file NORM_LOG2] [-calccov] [-hg38]  
                                   -a ANALYSIS [-tf TF]  
`

