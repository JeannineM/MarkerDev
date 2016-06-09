#!/usr/bin/env bash
#usage time ./gvcf_cal_Mrz-17.sh bam325.list > gvcf_cal_Mrz-17.$(date +%d-%m).log.txt 2> gvcf_cal_Mrz-17.$(date +%d-%m).err.txt
# START in sm11://data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/all331_Mrz-16/
# version March 2016
# realigned bam files

cd /data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/all331_Mrz-16/bamfiles/

OUTNAME=gvcffiles_Mrz-17
mkdir $OUTNAME

INPUT=$1

TARGETref=/data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/all336/refseq_gene236.fa
INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/interval_pic_300ampl_A.list
dbSNP=/data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/dbSNP_IDannot.vcf 


# echo "##################################################################################"
# echo '# Gatk: call all whole set of single-sample gVCFs'

# .. -ploidy
# -disableOptimizations \
# -mbq 20 \

# cat $INPUT | while read line; \
# do ID=$(sed 's/\..*$//' <<< $line); \
# 
# GenomeAnalysisTK.jar -T HaplotypeCaller \
# -R $TARGETref \
# -L $INTERV \
# -AR $INTERV \
# -D $dbSNP \
# -I ${ID}.vs236.coordSorted_RGh.realn.sorted.recal3.bam \
# -ERC GVCF \
# -dt NONE \
# -gt_mode DISCOVERY \
# -nda \
# -ploidy 2 \
# -mmq 44 \
# -mbq 20 \
# --max_alternate_alleles 2 \
# -bamWriterType CALLED_HAPLOTYPES -dontTrimActiveRegions -forceActive \
# -dontUseSoftClippedBases \
# -A Coverage -A MappingQualityZero -A FisherStrand -A QualByDepth -A DepthPerAlleleBySample \
# -A ChromosomeCounts -A VariantType -A DepthPerSampleHC -A AlleleCountBySample -A StrandBiasBySample -A GenotypeSummaries \
# -bamout ${OUTNAME}/${ID}.gatkHC.calledHaplotypes.bam \
# -o ${OUTNAME}/${ID}.gatkHC.raw.mmq44.mbq20.g.vcf; \
# 
# echo "DONE" $ID; \
# done
# echo "##################################################################################"

cd $OUTNAME

# ls *.gatkHC.raw.mmq44.mbq20.g.vcf > gvcf.list
#  
# echo 'gvcf.list: '; wc -l gvcf.list
# echo ''
# echo "###################################################################################"
# echo 'Merge gVCFs files -> get vcf file'
# echo "###################################################################################"

INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/amplicon_pic_300ampl.list
dbSNP=/data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/dbSNP_IDannot.vcf 
ID=sample325

# ls *.gatkHC.calledHaplotypes.bam | sort -n > bam.list
# time bamtools merge -list bam.list -out ${ID}.gatkHC.calledHaplotypes.bam
# samtools index ${ID}.gatkHC.calledHaplotypes.bam
# 
# time GenomeAnalysisTK.jar -T GenotypeGVCFs \
# -R $TARGETref \
# -L $INTERV \
# -D $dbSNP \
# -V gvcf.list \
# -dt NONE \
# -nt 8 \
# -stand_call_conf 30 \
# -stand_emit_conf 10 \
# -nda \
# --max_alternate_alleles 20 \
# -A Coverage -A MappingQualityZero -A FisherStrand -A QualByDepth -A DepthPerAlleleBySample -A InbreedingCoeff \
# -A ChromosomeCounts -A VariantType -A DepthPerSampleHC -A AlleleCountBySample -A StrandBiasBySample -A GenotypeSummaries \
# -o ${ID}.gatkgHC.raw.snps.indels_scc30_sec10_300ampl.vcf 
# 
# GenomeAnalysisTK.jar -T ValidateVariants \
# -R $TARGETref \
# -V ${ID}.gatkgHC.raw.snps.indels_scc30_sec10_300ampl.vcf
# 
# echo ''
# echo '###################################################################################'
# echo '# VQSR to GatK: Variant Recalibration'
# echo '###################################################################################'
# echo '# High Quality SNPs to be used in recalibration'
# 
# INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/interval_pic_184nuclear_good.list
# # = Processing 24464 bp from intervals 
# 
# GenomeAnalysisTK.jar -T VariantFiltration \
# -R $TARGETref \
# -dt NONE \
# -L $INTERV \
# -V ${ID}.gatkgHC.raw.snps.indels_scc30_sec10_300ampl.vcf \
# -filterName "HARD_TO_VALIDATE" -filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
# -filterName "HighFS" -filter "FS > 30.0" \
# -filterName "LowQD" -filter "QD < 2.0" \
# -filterName "LowCoverage" -filter "DP < 10" \
# -filterName "VeryLowQual" -filter "QUAL < 30.0" \
# -filterName "LoQual" -filter "QUAL > 30.0 && QUAL < 50.0" \
# -filterName "medQual" -filter "QUAL < 1000" \
# -filterName "noSNP" -filter "AF=='1.00'" \
# -cluster 3 -window 5 \
# -o ${ID}.gatkgHC.snps.indels_scc30_sec10_9filters_184nc.vcf
# 
# ### grep PASS only
# 
# cat ${ID}.gatkgHC.snps.indels_scc30_sec10_9filters_184nc.vcf | grep 'PASS\|^#' > highQualSNPs_184nc.vcf
# rm ${ID}.gatkgHC.snps.indels_scc30_sec10_9filters_184nc.vcf*
# 
# echo ''
# echo "###################################################################################"
# echo 'Variant Recalibration on SNPs gVCFcalls'
# echo "###################################################################################"
# 
# # -an InbreedingCoeff
# # --ignore_filter HARD_TO_VALIDATE \
# # --ignore_filter LowQD \
# # --ignore_filter LoQual \
# # --ignore_filter VeryLowQual \
# 
# INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/amplicon_pic_300ampl.list
# 
# GenomeAnalysisTK.jar -T VariantRecalibrator \
# -R $TARGETref \
# -dt NONE \
# -L $INTERV \
# -input ${ID}.gatkgHC.raw.snps.indels_scc30_sec10_300ampl.vcf \
# -resource:targetSet,known=true,training=true,truth=true,prior=10.0 /scratch/nine/cov_musagenesAug14/bangen/snpcalling_09-04/target.vcf \
# -resource:highQual,known=true,training=true,truth=true,prior=10.0 highQualSNPs_184nc.vcf \
# -resource:dbSNPs,known=true,training=false,truth=false,prior=2.0 /scratch/nine/cov_musagenesAug14/bangen/snpcalling_09-04/dbSNP_236genes.vcf \
# -an QD -an MQ -an MQRankSum -an ReadPosRankSum \
# -mode SNP \
# --maxGaussians 4 \
# --minNumBadVariants 5000 \
# -recalFile VQSR.recal \
# -tranchesFile VQSR.tranches \
# -rscriptFile VQSR.plots.R \
# -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0
# 
# GenomeAnalysisTK.jar -T ApplyRecalibration \
# -R $TARGETref \
# -dt NONE \
# -L $INTERV \
# -input ${ID}.gatkgHC.raw.snps.indels_scc30_sec10_300ampl.vcf \
# --ts_filter_level 90.0 \
# -mode SNP \
# -tranchesFile VQSR.tranches \
# -recalFile VQSR.recal \
# -o ${ID}.gatkgHC.snps.indels_scc30_sec10_300ampl.SNPrecal90.vcf
# 
# echo ''
# echo '###################################################################################'
# echo 'Hard filtering'
# echo "###################################################################################"
# 
# # undefined -filterName "LowQD" -filter "QD < 2.0" \
# INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/amplicon_pic_300ampl.list
# 
# GenomeAnalysisTK.jar -T VariantFiltration \
# -R $TARGETref \
# -L $INTERV \
# -dt NONE \
# -V ${ID}.gatkgHC.snps.indels_scc30_sec10_300ampl.SNPrecal90.vcf \
# -filterName "HARD_TO_VALIDATE" -filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
# -filterName "HighFS" -filter "FS > 30.0" \
# -filterName "LowCoverage" -filter "DP < 10" \
# -filterName "VeryLowQual" -filter "QUAL < 30.0" \
# -filterName "noSNP" -filter "AF=='1.00'" \
# -o ${ID}.gatkgHC.snps.indels_scc30_sec10_300ampl.SNPrecal90_7filt.vcf
# 
# echo '###################################################################################'
# echo '# Select variants'
# echo '# within the amplicons (length incl primers)'
# echo '# remove 4 poor samples with more than 12% no-call + 2 triploids = 325 samples'
# echo '# remove sites with more than 30% no-call (325*2 * (1-0.3))'
# echo '###################################################################################'
# 
# INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/interval_pic_300ampl_A.list

# -sn SAMPLE_A_PARC \
# -sn SAMPLE_B_ACTG \
# -xl_sf /data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/all336/gvcffiles/poor_sample.list \
# -ef exclude filtered

# -sn $samplelist \
# nocall=$1; echo "AN = $nocall"

nocall=455

INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/interval_pic_300ampl_A.list 

# GenomeAnalysisTK.jar -T SelectVariants \
# -R $TARGETref \
# -V ${ID}.gatkgHC.snps.indels_scc30_sec10_300ampl.SNPrecal90_7filt.vcf \
# -L $INTERV \
# -select "AN > $nocall" \
# -env \
# -dt NONE \
# -o ${ID}.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455.vcf

VCFinput=${ID}.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455.vcf

# GenomeAnalysisTK.jar -T VariantsToTable \
# -R $TARGETref \
# -L $INTERV \
# -V $VCFinput \
# -F CHROM -F POS -F ID -F REF -F ALT \
# -F NO-CALL -F HOM-REF -F HOM-VAR -F HET \
# -F AC -F AF -F AN -F VariantType \
# -GF GT \
# -o ${ID}.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455.tab

echo "###################################################################################"
echo '# Explore filters for the gvcf files' $VCFinput 
echo "###################################################################################"

echo 'Applied filters: '

VT=$(grep -c -v '^#' $VCFinput); echo '##Total number of called variants: ' $VT
VF=$(grep  -v '^#' $VCFinput | grep -vc "PASS"); echo '##Total number of filtered variants: ' $VF
VP=$(grep  -v '^#' $VCFinput | grep -c "PASS"); echo '##Total number of PASS variants: ' $VP
Genes=$(grep -v '^#' $VCFinput | cut -f 1| sort | uniq | wc -l); echo '##Occurring in  number of genes:' $Genes 

echo ''
echo 'Variant types:'
grep -v '^#' ${ID}.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455.tab | cut -f 13 | cut -f 1 -d "." | sort| uniq -c

echo ''
echo 'Known from transcriptome or new variants:'
grep -v '^#' $VCFinput | grep PASS | cut -f 3 | sort | uniq -c


echo ''
echo "###################################################################################"
echo '########## Final vcf file #######################################################'
echo "###################################################################################"
echo ''

cat ${ID}.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455.vcf | grep 'PASS\|^#' > ${ID}.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455_grepas.vcf

echo '-> ' ${ID}'.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455.vcf'
echo '-> ' ${ID}'.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455_grepas.vcf'
echo ''
echo "###################################################################################"
echo '# Merge phased variants to MNP and select S-M-NPs
echo "###################################################################################"
echo ''

TARGETref=/data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/all336/refseq_gene236.fa
INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/interval_pic_300ampl_A.list 
ID=sample325
nocall=455

GenomeAnalysisTK.jar -T ReadBackedPhasing \
-R $TARGETref \
-I sample325.gatkHC.calledHaplotypes.bam \
-V ${ID}.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455.vcf \
-L $INTERV \
-dt NONE \
--phaseQualityThresh 20.0 \
-cacheWindow 300 \
-maxDistMNP 200 \
-maxSites 300 \
-enableMergeToMNP \
-o ${ID}.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455.phasedMNP.vcf

echo ''
echo "###################################################################################"
echo '# Select all variants for subset of 71 samples and diversity measures
echo "###################################################################################"
echo ''

TARGETref=/data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/all336/refseq_gene236.fa
INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/interval_pic_300ampl_A.list 
ID=sample325
nocall=99


GenomeAnalysisTK.jar -T SelectVariants \
-R $TARGETref \
-V ${ID}.gatkgHC.snps.indels_scc30_sec10_300interv.SNPrecal90_7filt_AN455.phasedMNP.vcf \
-L $INTERV \
-dt NONE \
-sf markerPaperSample.list \
-env \
-trimAlternates \
-select "AN > $nocall" \
-o samples71.gatkgHC.snps.indels.scc30_sec10_300interv.SNPrecal90_7filt_AN99.phasedMNP.vcf

ID=samples71

GenomeAnalysisTK.jar -T VariantFiltration \
-R $TARGETref \
-L $INTERV \
-dt NONE \
-V ${ID}.gatkgHC.snps.indels.scc30_sec10_300interv.SNPrecal90_7filt_AN99.phasedMNP.vcf \
-filterName "noSNP" -filter "AF=='1.00'" \
-o ${ID}.gatkgHC.snps.indels.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.vcf


echo ''
echo "###################################################################################"
echo '# Select biSNPs for subset of 71 samples and data conversion
echo "###################################################################################"
echo ''

# remove sites with more than 30% no-call (71*2 * (1-0.3))'
nocall=99
# -select "AN > $nocall" \
# -selectType MNP

GenomeAnalysisTK.jar -T SelectVariants \
-R $TARGETref \
-L $INTERV \
-V ${ID}.gatkgHC.snps.indels.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.vcf \
-sf markerPaperSample.list \
-env -ef \
-trimAlternates \
-selectType SNP \
-select "AN > $nocall" \
-o samples71.gatkgHC.SNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.vcf

grep -v VariantType=MULTIALLELIC samples71.gatkgHC.SNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.vcf > samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.vcf

GenomeAnalysisTK.jar -T ValidateVariants \
-R $TARGETref \
-V samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.vcf

VCFinput=samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.vcf

/data/home/btw680/software/vcftools_0.1.12b/bin/vcftools --vcf $VCFinput --012 --out samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP

GenomeAnalysisTK.jar -T VariantsToTable \
-R $TARGETref \
-L $INTERV \
-V $VCFinput \
-F CHROM -F POS -F ID -F REF -F ALT \
-F NO-CALL -F HOM-REF -F HOM-VAR -F HET \
-F AC -F AF -F AN -F VariantType \
-GF GT \
-o samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.tab

echo ''
echo "###################################################################################"
echo '# Select biSNPs for subset of 71 samples GOOD GENES
echo "###################################################################################"
echo ''

INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/interval_pic_198_JMgood.list

VCFinput=samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.vcf

GenomeAnalysisTK.jar -T SelectVariants \
-R $TARGETref \
-V $VCFinput \
-L $INTERV \
-env -ef \
-trimAlternates \
-o samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L198good.vcf

GenomeAnalysisTK.jar -T ValidateVariants \
-R $TARGETref \
-V samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L198good.vcf

VCFinput=samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L198good.vcf

/data/home/btw680/software/vcftools_0.1.12b/bin/vcftools --vcf $VCFinput --012 --out samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L198good

GenomeAnalysisTK.jar -T VariantsToTable \
-R $TARGETref \
-L $INTERV \
-V $VCFinput \
-F CHROM -F POS -F ID -F REF -F ALT \
-F NO-CALL -F HOM-REF -F HOM-VAR -F HET \
-F AC -F AF -F AN -F VariantType \
-GF GT \
-o samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L198good.tab

echo ''
echo "###################################################################################"
echo '# Select biSNPs for subset of 71 samples RBGE GENES: 207
echo "###################################################################################"
echo ''

INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/amplicon_pic_207sub_RBGE.list

VCFinput=samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.vcf

GenomeAnalysisTK.jar -T SelectVariants \
-R $TARGETref \
-V $VCFinput \
-L $INTERV \
-env -ef \
-trimAlternates \
-o samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L207RBGE.vcf

GenomeAnalysisTK.jar -T ValidateVariants \
-R $TARGETref \
-V samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L207RBGE.vcf

VCFinput=samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L207RBGE.vcf

/data/home/btw680/software/vcftools_0.1.12b/bin/vcftools --vcf $VCFinput --012 --out samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L207RBGE

GenomeAnalysisTK.jar -T VariantsToTable \
-R $TARGETref \
-L $INTERV \
-V $VCFinput \
-F CHROM -F POS -F ID -F REF -F ALT \
-F NO-CALL -F HOM-REF -F HOM-VAR -F HET \
-F AC -F AF -F AN -F VariantType \
-GF GT \
-o samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L207RBGE.tab

echo ''
echo "###################################################################################"
echo '# FastStructure Input vcf ntDNA only
echo "###################################################################################"
echo ''

cd /data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/all331_Mrz-16/bamfiles/gvcffiles_Mrz-17/

TARGETref=/data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/all336/refseq_gene236.fa
INTERV=/data/scratch/btw680/analyses/reseq_ana/all_files/interval_pic_285nt_genes.list
VCFinput=samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.vcf

GenomeAnalysisTK.jar -T SelectVariants \
-R $TARGETref \
-V $VCFinput \
-L $INTERV \
-env -ef \
-trimAlternates \
-o samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L285nt_genes.vcf

GenomeAnalysisTK.jar -T ValidateVariants \
-R $TARGETref \
-V samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L285nt_genes.vcf

VCFinput=samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L285nt_genes.vcf
/data/home/btw680/software/vcftools_0.1.12b/bin/vcftools --vcf $VCFinput --012 --out samples71.gatkgHC.biSNPs.scc30_sec10_300interv.SNPrecal90_8filt_AN99.phasedMNP.L285nt_genes

echo ''
echo "###################################################################################"
echo "DONE:" $(date +%d-%m) $(date +%H:%M)
echo "###################################################################################"
echo ''
