#! /bin/bash

# Calculate Accessibility Scores
##########################################
# 1 accesibility by high-frequency signal amplitude, corrected by number of binding sites
#Activity score calculations
# (positional) command args:
#  1) Sample Folder
#  2) Sample name
echo "Calculate HighRange Frequency Amplitude"
mkdir ${1}/AccessibilityOutput
for i in ${1}/TranscriptionFactors/GTRD_ChIP_Only/*.tss
do
    tf=`basename $i`
    Rscript ./Scripts/activity_score_highfreq_signal_range.R $i $tf ${1}/AccessibilityOutput/${2}_AccessibilityAllHighFreq.txt
done

echo "Adjust amplitudes by loess smoothing"
Rscript ./Scripts/loess_adjust_high_freq.R ${1}/AccessibilityOutput/${2}_AccessibilityAllHighFreq.txt ${1}/AccessibilityOutput/${2}_AccessibilityAllHighFreqAdjusted_tmp.txt ${1}/AccessibilityOutput/${2}_loess_highfreq_all.png


cat ${1}/AccessibilityOutput/${2}_AccessibilityAllHighFreqAdjusted_tmp.txt >> ${1}/AccessibilityOutput/${2}_AccessibilityAllHighFreqAdjusted.txt
rm ${1}/AccessibilityOutput/${2}_AccessibilityAllHighFreqAdjusted_tmp.txt


##########################################
# 2 accessibility by signal amplitude of signal after deconstruction and detrending from wavelet analysis
# scoring pipeline using wavelets
mkdir ${1}/AccessibilityOutput/WaveletOut
echo -e "TranscriptionFactor\tBindingSites\tMaxPower_Period\tMaxPower\tMaxPower_PeriodInRange\tMaxPower_inRange\tReconstructAmplitude\tpower_reconstruct\tsum_ampl_reconstruct" > ./${2}/${2}_AccessibilityAllWavelets.txt
for i in ${1}/TranscriptionFactors/GTRD_ChIP_Only/*.tss
do
   tf=`basename $i`
   echo $tf
   Rscript ./Scripts/activityScore_wavelet.R $i $tf ${1}/AccessibilityOutput/WaveletOut/$tf ${1}/AccessibilityOutput/${2}_AccessibilityAllWavelets.txt
done

Rscript ./Scripts/loess_adjust_wavelets.R ${1}/AccessibilityOutput/${2}_AccessibilityAllWavelets.txt ${1}/AccessibilityOutput/${2}_AccessibilityAllWaveletsAdjusted.txt ${1}/AccessibilityOutput/${2}_loess_wavelet_all.png


##########################################
# 3 accesibility by high-frequency signal amplitude of 1000 binding sites of each tF (that has > 1000 binding sites defined)
#     this signal is not corrected by LOESS since the noise should be equal as there are an equal amount of binding sites analyzed
#Activity score calculations
# (positional) command args:
#  1) Folder with GTRD TSS files
#  2) Sample name
echo "Calculate HighRange Frequency Amplitude"
for i in ${1}/TranscriptionFactors/GTRD_ChIP_Only_1000sites/*.tss
do
    tf=`basename $i`
    Rscript ./Scripts/activity_score_highfreq_signal_range.R $i $tf ${1}/AccessibilityOutput/${2}_Accessibility1KSites.txt
done

#"No loess smoothing is done, but ranks are calculated"
Rscript ./Scripts/ranks_1K.R ${1}/AccessibilityOutput/${2}_Accessibility1KSites.txt ${1}/AccessibilityOutput/${2}_Accessibility1KSitesAdjusted_tmp.txt

echo -e "TF_profile\tTFBindingSites\tHighFreqRange\tPeaks\tMeanPeakDistance\tMedianPeakDistance\tHighFreqRange_rank" > ${1}/AccessibilityOutput/${2}_Accessibility1KSitesAdjusted.txt
cat ${1}/AccessibilityOutput/${2}_Accessibility1KSitesAdjusted_tmp.txt >> ${1}/AccessibilityOutput/${2}_Accessibility1KSitesAdjusted.txt
rm ${1}/AccessibilityOutput/${2}_Accessibility1KSitesAdjusted_tmp.txt


