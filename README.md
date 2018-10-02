# TranscriptionFactorProfiling
Profiling of transcription factor binding sites in cell-free DNA

Raw sequencing data is trimmed to contain bases 53-113 (i.e. the central 60bp in a hypothetical 166bp fragment), PCR duplicates are removed and trimmed read is aligned.

In a second step, coverage data around transcription factor binding sites are averaged and the activity for each transcription factor is calculated.

Transcription factor binding sites are taken from GTRD database

Software needed:  
* bedtools  
* bwa  
* samtools  
