#!/usr/bin/env bash
#usage ./Snp_calling_MusA.sh > screen_$(date +%d-%m).log.txt 2> screen_$(date +%d-%m).err.txt
# 
#input files are the (realigned, recalibrated, sorted and duplicate reads removed) bam 
#files from realignSWA{1,2,3}.sh

cd /scratch/nine/cov_musagenesAug14/bangen
mkdir snpcalling_$(date +%d-%m); cd snpcalling_$(date +%d-%m)

TARGETref=refseqs_1046_manedt.fa
OUTNAME=BBSWA1-3vs1046MusA

INNAME1=SWA1vs1046genes
INNAME2=SWA2vs1046genes
INNAME3=SWA3vs1046genes

cp ../${INNAME1}/${INNAME1}.coordSorted_RGh.raln.recal3.sorted.Srmdup.ba* ./
cp ../${INNAME2}/${INNAME2}.coordSorted_RGh.raln.recal3.sorted.Srmdup.ba* ./
cp ../${INNAME3}/${INNAME3}.coordSorted_RGh.raln.recal3.sorted.Srmdup.ba* ./
cp /scratch/nine/cov_musagenesAug14/RG.txt ./

cp /scratch/nine/cov_musagenesAug14/$TARGETref ./

CreateSequenceDictionary.jar R=$TARGETref O=${TARGETref}.dict
cp ${TARGETref}.dict refseqs_1046_manedt.dict

samtools faidx $TARGETref

echo '###################################################################################'
echo '# Merge bam files'

samtools merge ${OUTNAME}.bammerged.bam -h RG.txt ${INNAME1}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam ${INNAME2}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam ${INNAME3}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam 
samtools index ${OUTNAME}.bammerged.bam 

cd /scratch/nine/cov_musagenesAug14/bangen/snpcalling_$(date +%d-%m)

GenomeAnalysisTK.jar -T RealignerTargetCreator -R $TARGETref -I ${OUTNAME}.bammerged.bam -o ${OUTNAME}.bammerged.intervals
GenomeAnalysisTK.jar -T IndelRealigner -R $TARGETref -I ${OUTNAME}.bammerged.bam -targetIntervals ${OUTNAME}.bammerged.intervals -o ${OUTNAME}.bammerged.realn.bam
samtools index ${OUTNAME}.bammerged.realn.bam

samtools sort ${OUTNAME}.bammerged.realn.bam ${OUTNAME}.bammerged.realn.sorted
samtools index ${OUTNAME}.bammerged.realn.sorted.bam

echo "##################################################################################"
echo "Get coverage for merged bam file (all 3 species)"

bamToBed -i ${OUTNAME}.bammerged.realn.sorted.bam > ${OUTNAME}.bammerged.realn.sorted.bed
samtools depth -b ${OUTNAME}.bammerged.realn.sorted.bed ${OUTNAME}.bammerged.realn.sorted.bam > ${OUTNAME}.bammerged.realn.sorted_samdepth.txt
Rscript /scratch/nine/cov_musagenesAug14/coverage.R ${OUTNAME}.bammerged.realn.sorted_samdepth.txt ${OUTNAME}.bammerged.realn.sorted_cov.txt

echo '###################################################################################'
echo '# Gatk: call all low confidence SNPs'

### Haplotype caller 

GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $TARGETref \
-I ${INNAME1}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
-I ${INNAME2}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
-I ${INNAME3}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
-gt_mode DISCOVERY \
-stand_call_conf 3 \
-stand_emit_conf 4 \
-minPruning 5 \
-nda \
-ploidy 2 \
-A Coverage -A MappingQualityZero -A FisherStrand -A QualByDepth -A DepthPerAlleleBySample \
-A ChromosomeCounts -A VariantType \
-o ${OUTNAME}.gatkHC_LQ.raw.snps.indels_scc3_sec4.vcf &

pidHC1=$!
echo $pidHC1
echo ''
echo '###################################################################################'
echo '# Gatk: call high confidence SNPs for calibration'

### Haplotype caller

# GenomeAnalysisTK.jar -T HaplotypeCaller \
# -R $TARGETref \
# -I ${INNAME1}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
# -I ${INNAME2}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
# -I ${INNAME3}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
# -gt_mode DISCOVERY \
# -stand_call_conf 20 \
# -stand_emit_conf 20 \
# -minPruning 5 \
# -nda \
# -ploidy 2 \
# -o ${OUTNAME}.gatkHC_HQ.raw.snps.indels_scc20_sec20.vcf &

GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $TARGETref \
-I ${INNAME1}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
-I ${INNAME2}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
-I ${INNAME3}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
-gt_mode DISCOVERY \
-stand_call_conf 20 \
-stand_emit_conf 20 \
-minPruning 5 \
-nda \
-ploidy 2 \
-A Coverage -A MappingQualityZero -A FisherStrand -A QualByDepth -A DepthPerAlleleBySample \
-A ChromosomeCounts -A VariantType \
-o ${OUTNAME}.gatkHC_LQ.raw.snps.indels_scc20_sec20.vcf &



pidHC2=$!
echo 'pidHC2: ' $pidHC2
echo ''
echo '###################################################################################'
# samtools variant calling
# 
# echo 'samtools multi-sample raw snp calling'
# 
# echo 'options samtools mpileup: \
# - no Indel calling \
# - raw read depth coverage per sample \
# - generate (uncompressed) bcf (=genotype likelihood) output \
# - extended BAQ for higher sensitivity \
# - C50 \
# - maximum per BAM depth 1000 '
# 
# samtools mpileup -DEugI -C50 -d1000 -f $TARGETref \
# ${INNAME1}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
# ${INNAME2}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam \
# ${INNAME3}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam | /share/apps/sbcs/bcftools/0.1.17/bin/bcftools view -bvcg - > ${OUTNAME}.sam.raw_C50_E.bcf
# 
# echo 'applied filters for output sam vcf: \
# - print filtered variants \
# - min read depth 5 '
# 
# /share/apps/sbcs/bcftools/0.1.17/bin/bcftools view ${OUTNAME}.sam.raw_C50_E.bcf | /share/apps/sbcs/bcftools/0.1.17/bin/vcfutils.pl varFilter -p -d5 > ${OUTNAME}.sam.raw.C50_E_d5.vcf
# 
# ### Select SNPs only
# 
# GenomeAnalysisTK.jar -T SelectVariants \
# -R $TARGETref \
# -V ${OUTNAME}.sam.raw.C50_E_d5.vcf \
# -selectType SNP \
# -o ${OUTNAME}.sam.snps.C50_E_d5.vcf 
# 
# ### hard filtering
# 
# echo ' SAM hard filters '
# 
# GenomeAnalysisTK.jar -T VariantAnnotator \
# -R $TARGETref \
# -V ${OUTNAME}.sam.snps.C50_E_d5.vcf \
# -I ${OUTNAME}.bammerged.realn.sorted.bam \
# -A Coverage -A MappingQualityZero -A FisherStrand -A QualByDepth -A DepthPerAlleleBySample \
# -A ChromosomeCounts -A VariantType \
# -o ${OUTNAME}.sam.snps.C50_E_d5.anno.vcf 
# 
# GenomeAnalysisTK.jar -T VariantFiltration \
# -R $TARGETref \
# --invalidatePreviousFilters \
# -V ${OUTNAME}.sam.snps.C50_E_d5.anno.vcf \
# -filterName "HARD_TO_VALIDATE" -filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
# -filterName "HighFS" -filter "FS > 30.0" \
# -filterName "LowQD" -filter "QD < 2.0" \
# -filterName "LowCoverage" -filter "DP < 5" \
# -filterName "VeryLowQual" -filter "QUAL < 30.0" \
# -filterName "LoQual" -filter "QUAL > 30.0 && QUAL < 50.0" \
# -filterName "noSNP" -filter "AF=='1.00'" \
# -o ${OUTNAME}.sam.snps.C50_E_d5.anno_7filters.vcf
#    
# echo ' '
# echo 'DONE samSNPs'
# 
echo '###################################################################################'

wait $pidHC1
wait $pidHC2

echo '###################################################################################'
echo 'Gatk: LowQuality Select SNPs only'

GenomeAnalysisTK.jar -T SelectVariants \
-R $TARGETref \
-V ${OUTNAME}.gatkHC_LQ.raw.snps.indels_scc3_sec4.vcf \
-selectType SNP \
-o ${OUTNAME}.gatkHC_LQ.raw.snps_scc3_sec4.vcf 

### hard filtering

GenomeAnalysisTK.jar -T VariantFiltration \
-R $TARGETref \
-V ${OUTNAME}.gatkHC_LQ.raw.snps_scc3_sec4.vcf \
-filterName "HARD_TO_VALIDATE" -filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
-filterName "HighFS" -filter "FS > 30.0" \
-filterName "LowQD" -filter "QD < 2.0" \
-filterName "LowCoverage" -filter "DP < 5" \
-filterName "VeryLowQual" -filter "QUAL < 30.0" \
-filterName "LoQual" -filter "QUAL > 30.0 && QUAL < 50.0" \
-filterName "noSNP" -filter "AF=='1.00'" \
-o ${OUTNAME}.gatkHC_LQ.snps_scc3_sec4_7filters.vcf


### Gatk: HighQuality Select SNPs only

GenomeAnalysisTK.jar -T SelectVariants \
-R $TARGETref \
-V ${OUTNAME}.gatkHC_HQ.raw.snps.indels_scc20_sec20.vcf \
-selectType SNP \
-o ${OUTNAME}.gatkHC_HQ.raw.snps_scc20_sec20.vcf 

### hard filtering

GenomeAnalysisTK.jar -T VariantFiltration \
-R $TARGETref \
-V ${OUTNAME}.gatkHC_HQ.raw.snps_scc20_sec20.vcf \
-window 35 -cluster 3 \
-filterName "HARD_TO_VALIDATE" -filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
-filterName "HighFS" -filter "FS > 30.0" \
-filterName "LowQD" -filter "QD < 2.0" \
-filterName "LowCoverage" -filter "DP < 10" \
-filterName "VeryLowQual" -filter "QUAL < 30.0" \
-filterName "LoQual" -filter "QUAL > 30.0 && QUAL < 50.0" \
-filterName "medQual" -filter "QUAL < 500" \
-filterName "NAllsamples" -filter "AN < '6'" \
-filterName "noSNP" -filter "AF=='1.00'" \
-o ${OUTNAME}.gatkHC_HQ.snps_scc20_sec20_10filters.vcf


### grep PASS only

cat ${OUTNAME}.gatkHC_HQ.snps_scc20_sec20_10filters.vcf | grep 'PASS\|^#' > highQualSNPS.vcf

echo '###################################################################################'
# #5 Combine & compare Variants
# 
# GenomeAnalysisTK.jar -T CombineVariants \
# -R $TARGETref \
# -V:Gatk ${OUTNAME}.gatkHC_LQ.snps_scc3_sec4_7filters.vcf \
# -V:Sam  ${OUTNAME}.sam.snps.C50_E_d5.anno_7filters.vcf \
# -minN 1 \
# -env \
# -genotypeMergeOptions PRIORITIZE \
# -priority Gatk,Sam \
# -o ${OUTNAME}.GatkSam.merged.SNPs.vcf
# 
# GenomeAnalysisTK.jar -T VariantAnnotator \
# -R $TARGETref \
# -V ${OUTNAME}.GatkSam.merged.SNPs.vcf \
# -I ${OUTNAME}.bammerged.realn.sorted.bam \
# -A Coverage -A MappingQualityZero -A FisherStrand -A QualByDepth -A DepthPerAlleleBySample \
# -A VariantType -A MappingQualityRankSumTest -A ReadPosRankSumTest -A ChromosomeCounts \
# -o ${OUTNAME}.GatkSam.merged.SNPs.anno.vcf
# 
# GenomeAnalysisTK.jar -T VariantFiltration \
# -R $TARGETref \
# -V ${OUTNAME}.GatkSam.merged.SNPs.anno.vcf \
# --invalidatePreviousFilters \
# -filterName "HARD_TO_VALIDATE" -filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
# -filterName "HighFS" -filter "FS > 30.0" \
# -filterName "LowQD" -filter "QD < 2.0" \
# -filterName "LowCoverage" -filter "DP < 10" \
# -filterName "VeryLowQual" -filter "QUAL < 30.0" \
# -filterName "LoQual" -filter "QUAL > 30.0 && QUAL < 50.0" \
# -filterName "noSNP" -filter "AF=='1.00'" \
# -o ${OUTNAME}.GatkSam.merged.SNPs.anno_7filters.vcf
# 

echo '###################################################################################'
echo '# VQSR to GatK'
### VariantRecalibration
INPUT=BBSWA1-3vs1046MusA.gatkHC_LQ.snps_scc3_sec4_7filters.vcf

GenomeAnalysisTK.jar -T VariantRecalibrator \
-R $TARGETref \
-input $INPUT \
-resource:concordantSet,known=true,training=true,truth=true,prior=10.0 highQualSNPS.vcf \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
-mode SNP \
--maxGaussians 4 \
-recalFile VQSR.recal \
-tranchesFile VQSR.tranches \
-rscriptFile VQSR.plots.R \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
--ignore_filter HARD_TO_VALIDATE \
--ignore_filter LoQual \
--ignore_filter SNPcluster \
--ignore_filter medQual

# --ignore_filter HighFS \
# --ignore_filter LowQD \
# --ignore_filter LowCoverage \
# --ignore_filter VeryLowQual \


### Apply Recalibration

GenomeAnalysisTK.jar -T ApplyRecalibration \
-R $TARGETref \
-input BBSWA1-3vs1046MusA.gatkHC_LQ.snps_scc3_sec4_7filters.vcf \
--ts_filter_level 99.0 \
-mode SNP \
--ignore_filter HARD_TO_VALIDATE \
--ignore_filter LoQual \
--ignore_filter SNPcluster \
--ignore_filter medQual \
-tranchesFile VQSR.tranches \
-recalFile VQSR.recal \
-o BBSWA1-3vs1046MusA.gatkHC_LQ.snps_scc3_sec4_7filters.recal.vcf

echo '###################################################################################'
echo '# Haplotype Phasing with high confidence SNPs'

### ReadBackedPhasing

GenomeAnalysisTK.jar -T ReadBackedPhasing \
-R $TARGETref \
-I ${OUTNAME}.bammerged.realn.sorted.bam \
-V ${OUTNAME}.gatkHC_LQ.snps_scc3_sec4_7filters.recal.vcf \
--phaseQualityThresh 20.0 \
-o ${OUTNAME}.gatkHC_LQ.snps_scc3_sec4_7filters.recal.phased.vcf

### CalculateGenotypePosterior ## only valid for over 10 samples

GenomeAnalysisTK.jar -T CalculateGenotypePosteriors \
-R $TARGETref \
-V ${OUTNAME}.gatkHC_LQ.snps_scc3_sec4_7filters.recal.phased.vcf \
-useAC \
-o ${OUTNAME}.gatkHC_LQ.snps_scc3_sec4_7filters.recal.phased.withPost.AC.vcf


echo ''
echo '###########  Final vcf file ########################################################################'
echo ''

# grep PASS

cat ${OUTNAME}.gatkHC_LQ.snps_scc3_sec4_7filters.recal.phased.withPost.AC.vcf | grep 'PASS\|^#' > ${OUTNAME}.gatkHC_LQ.snps_scc3_sec4_7filters.recal.phased.withPost.AC_grepas.vcf

echo 'Outfiles: ' 
echo '-> ' ${OUTNAME}'.gatkHC_LQ.snps_scc3_sec4_7filters.recal.phased.withPost.AC.vcf'
echo '-> ' ${OUTNAME}'.gatkHC_LQ.snps_scc3_sec4_7filters.recal.phased.withPost.AC_grepas.vcf'




