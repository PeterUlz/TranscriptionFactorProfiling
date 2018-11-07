# TranscriptionFactorProfiling
Profiling of transcription factor binding sites in cell-free DNA

Raw sequencing data is trimmed to contain bases 53-113 (i.e. the central 60bp in a hypothetical 166bp fragment), PCR duplicates are removed and trimmed read is aligned.

In a second step, coverage data around transcription factor binding sites are averaged and the activity for each transcription factor is calculated.

Transcription factor binding sites are taken from GTRD database (http://gtrd.biouml.org/)

## Preparation:
Software needed:
* bedtools (v2.24.0)
* bwa (0.7.4-r385)
* samtools (v 1.3)
* fastx_trimmer (0.0.13)
* picard (as jar file)

Binaries should be put in ./Software/

## Trim and align
FastQ files are trimmed stripping away anything between base 53-113. This should still contain the most central 60bp of a hypothetical 166 bp fragment and thus be the
portion of the fragment that is closest to the nucleosome dyad. 

`
usage: trim_and_align.py [-h] -fq FASTQ_FILE -s NAME [-o OUTDIR] [-k]
                         [-t THREADS] [-step START_STEP]
`

## Trim from BAM file
Single-end aligned BAM files can also be diretly manipulated uing trim_from_bam_single_end.py
This is useful when the FastQ files are short than 113 bp.
However, in order to increase speed the BAM loses a lot of details: CIGAR strings, sequence is replaced with all "A"s (in order not to
look up the reference sequence)
The trimmed BAM is still good for coverage analyses, though.  
`
usage: trim_from_bam_single_end.py [-h] -i IN_FILE -o OUT_FILE -l LENGTH
`

## Calculate copy-number alterations
This is used in order to correct for local coverage differences due to copy-number alterations. 
The output of this is directly fed into the subsequent program, where coverage 
One option to do this is Plasma-Seq (https://github.com/PeterUlz/PlasmaSeq).
Other options could include: ichorCNA (https://github.com/broadinstitute/ichorCNA) or qDNAseq (https://bioconductor.org/packages/release/bioc/html/QDNAseq.html)

However, make sure the input file for the following steps looks like this:  
`
<chr>TAB<start>TAB<end>TAB<Log2-ratio>
`  
This corresponds to the .segments data produced by Plasma-Seq


## Run coverage analysis over all transcription factors in GTRD
In the next step, average coverages around defined transcription factor binding sites are calculated  

`
usage: run_tf_analyses_from_bam.py [-h] -b BAM_FILE -o NAME  
                                   [-cov MEAN_COVERAGE] [-ylimit YLIMIT]  
                                   [-norm-file NORM_LOG2] [-calccov] [-hg38]  
                                   -a ANALYSIS [-tf TF]  
`

## Calculate Accessibility per TF
After average coverages around TF binding sites have been calculated, we can estimate the accessibility by specifying the output folder of the step above and
the respective sample name  
`
./scoring_pipeline.sh <output_folder> <sample name>
`  



