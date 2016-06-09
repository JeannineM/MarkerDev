#!/usr/bin/env bash
#usage ./Snp_calling_Plastid.sh > screen_$(date +%d-%m).log.txt 2> screen_$(date +%d-%m).err.txt
# 
#input files are the (realigned, recalibrated, sorted and duplicated reads removed) bam 
#files from realignSWA{1,2,3}.sh

cd /scratch/nine/cov_musagenesAug14/plastid
mkdir snpcalling_$(date +%d-%m); cd snpcalling_$(date +%d-%m)

TARGETref=refseq84cp40mtPhoDcds.fa
OUTNAME=BBSWA1-3vs124pl

INNAME1=SWA1vs124plgenes
INNAME2=SWA2vs124plgenes
INNAME3=SWA3vs124plgenes

cp ../${INNAME1}/${INNAME1}.coordSorted_RGh.raln.recal3.sorted.Srmdup.ba* ./
cp ../${INNAME2}/${INNAME2}.coordSorted_RGh.raln.recal3.sorted.Srmdup.ba* ./
cp ../${INNAME3}/${INNAME3}.coordSorted_RGh.raln.recal3.sorted.Srmdup.ba* ./
cp /scratch/nine/cov_musagenesAug14/RG.txt ./

cp /scratch/nine/cov_musagenesAug14/$TARGETref ./

CreateSequenceDictionary.jar R=$TARGETref O=${TARGETref}.dict
cp ${TARGETref}.dict refseq84cp40mtPhoDcds.dict

samtools faidx $TARGETref

echo '###################################################################################'
echo '# merge bam files'

samtools merge ${OUTNAME}.bammerged.bam -h RG.txt ${INNAME1}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam ${INNAME2}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam ${INNAME3}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam 
samtools index ${OUTNAME}.bammerged.bam 

cd /scratch/nine/cov_musagenesAug14/plastid/snpcalling_$(date +%d-%m)

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

a=$(grep "NC_016740" ${OUTNAME}.bammerged.realn.sorted_samdepth.txt | wc -l)
echo 'length covered mt bases: ' $a
a=$(grep "NC_013991" ${OUTNAME}.bammerged.realn.sorted_samdepth.txt | wc -l)
	echo 'length covered cp bases: ' $a

echo ''
echo '###################################################################################'
echo '############## DIPLOID calling ####################################################'
echo '###################################################################################'
echo 'Gatk: call diploid SNPs'

### Haplotype caller 

GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $TARGETref \
-I ${OUTNAME}.bammerged.realn.sorted.bam \
-gt_mode DISCOVERY \
-dontUseSoftClippedBases \
-stand_call_conf 20 \
-stand_emit_conf 20 \
-minPruning 5 \
-nda \
-A Coverage -A MappingQualityZero -A FisherStrand -A QualByDepth -A DepthPerAlleleBySample \
-A ChromosomeCounts -A VariantType \
-ploidy 2 \
-o ${OUTNAME}.gatkHC.raw.snps.indels_scc20_sec20.vcf

echo '###################################################################################'
echo 'Gatk: Select SNPs only'

GenomeAnalysisTK.jar -T SelectVariants \
-R $TARGETref \
-V ${OUTNAME}.gatkHC.raw.snps.indels_scc20_sec20.vcf \
-selectType SNP \
-o ${OUTNAME}.gatkHC.raw.snps_scc20_sec20.vcf 

### hard filtering

GenomeAnalysisTK.jar -T VariantFiltration \
-R $TARGETref \
-V ${OUTNAME}.gatkHC.raw.snps_scc20_sec20.vcf  \
-filterName "HARD_TO_VALIDATE" -filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
-filterName "HighFS" -filter "FS > 30.0" \
-filterName "LowQD" -filter "QD < 2.0" \
-filterName "LowCoverage" -filter "DP < 5" \
-filterName "VeryLowQual" -filter "QUAL < 30.0" \
-filterName "LoQual" -filter "QUAL > 30.0 && QUAL < 50.0" \
-filterName "noSNP" -filter "AF=='1.00'" \
-o ${OUTNAME}.gatkHC.snps.scc20_sec20_7filters.vcf

echo ''
echo '###################################################################################'
echo 'Haplotype Phasing with high confidence SNPs'

### ReadBackedPhasing

GenomeAnalysisTK.jar -T ReadBackedPhasing \
-R $TARGETref \
-I ${OUTNAME}.bammerged.realn.sorted.bam \
-V ${OUTNAME}.gatkHC.snps.scc20_sec20_7filters.vcf \
--phaseQualityThresh 20.0 \
-o ${OUTNAME}.gatkHC.snps.scc20_sec20_7filters.phased.vcf

echo ''
echo '###################################################################################'
echo 'Split cp from mt into separate files'

# NC_013991.2 84 = cp
# NC_016740.1 40 = mt

cat ${OUTNAME}.gatkHC.snps.scc20_sec20_7filters.phased.vcf | grep -v 'NC_016740' > ${OUTNAME}.cp.2n.finalSnps.vcf
cat ${OUTNAME}.gatkHC.snps.scc20_sec20_7filters.phased.vcf | grep -v 'NC_013991' > ${OUTNAME}.mt.2n.finalSnps.vcf

cat ${OUTNAME}.cp.2n.finalSnps.vcf | grep "PASS\|^#" > ${OUTNAME}.cp.2n.finalSnps_grepas.vcf
cat ${OUTNAME}.mt.2n.finalSnps.vcf | grep "PASS\|^#" > ${OUTNAME}.mt.2n.finalSnps_grepas.vcf

GenomeAnalysisTK.jar -T ValidateVariants \
-R $TARGETref \
--variant ${OUTNAME}.cp.2n.finalSnps.vcf

GenomeAnalysisTK.jar -T ValidateVariants \
-R $TARGETref \
--variant ${OUTNAME}.mt.2n.finalSnps.vcf

echo '->' ${OUTNAME}'.cp.finalSnps.vcf '
VCFGREPFILEcp=${OUTNAME}.cp.2n.finalSnps_grepas.vcf
VCFFILEcp=${OUTNAME}.cp.2n.finalSnps.vcf

echo '->' ${OUTNAME}'.mt.finalSnps.vcf '
VCFGREPFILEmt=${OUTNAME}.mt.2n.finalSnps_grepas.vcf
VCFFILEmt=${OUTNAME}.mt.2n.finalSnps.vcf

echo '###################################################################################'
echo '############## HAPLOID calling ####################################################'
echo '###################################################################################'
echo 'Gatk: call haploid SNPs'

### Haplotype caller 

GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $TARGETref \
-I ${OUTNAME}.bammerged.realn.sorted.bam \
-gt_mode DISCOVERY \
-dontUseSoftClippedBases \
-stand_call_conf 20 \
-stand_emit_conf 20 \
-minPruning 5 \
-nda \
-A Coverage -A MappingQualityZero -A FisherStrand -A QualByDepth -A DepthPerAlleleBySample \
-A ChromosomeCounts -A VariantType \
-ploidy 1 \
-o ${OUTNAME}.gatkHC.raw.1n.snps.indels_scc20_sec20.vcf

echo '###################################################################################'
echo 'Gatk: Select haploid SNPs only'

GenomeAnalysisTK.jar -T SelectVariants \
-R $TARGETref \
-V ${OUTNAME}.gatkHC.raw.1n.snps.indels_scc20_sec20.vcf \
-selectType SNP \
-o ${OUTNAME}.gatkHC.raw.1n.snps_scc20_sec20.vcf 

### hard filtering

GenomeAnalysisTK.jar -T VariantFiltration \
-R $TARGETref \
-V ${OUTNAME}.gatkHC.raw.1n.snps_scc20_sec20.vcf  \
-filterName "HARD_TO_VALIDATE" -filter "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
-filterName "HighFS" -filter "FS > 30.0" \
-filterName "LowQD" -filter "QD < 2.0" \
-filterName "LowCoverage" -filter "DP < 5" \
-filterName "VeryLowQual" -filter "QUAL < 30.0" \
-filterName "LoQual" -filter "QUAL > 30.0 && QUAL < 50.0" \
-filterName "noSNP" -filter "AF=='1.00'" \
-o ${OUTNAME}.gatkHC.1n.snps.scc20_sec20_7filters.vcf

echo '###################################################################################'
echo 'Haplotype Phasing with high confidence SNPs'

### ReadBackedPhasing

GenomeAnalysisTK.jar -T ReadBackedPhasing \
-R $TARGETref \
-I ${OUTNAME}.bammerged.realn.sorted.bam \
-V ${OUTNAME}.gatkHC.1n.snps.scc20_sec20_7filters.vcf \
--phaseQualityThresh 20.0 \
-o ${OUTNAME}.gatkHC.1n.snps.scc20_sec20_7filters.phased.vcf

echo '###################################################################################'
echo 'Split cp from mt into separate files'

# NC_013991.2 84 = cp
# NC_016740.1 40 = mt

cat ${OUTNAME}.gatkHC.1n.snps.scc20_sec20_7filters.phased.vcf | grep -v 'NC_016740' > ${OUTNAME}.cp.1n.finalSnps.vcf
cat ${OUTNAME}.gatkHC.1n.snps.scc20_sec20_7filters.phased.vcf | grep -v 'NC_013991' > ${OUTNAME}.mt.1n.finalSnps.vcf

cat ${OUTNAME}.cp.1n.finalSnps.vcf | grep "PASS\|^#" > ${OUTNAME}.cp.1n.finalSnps_grepas.vcf
cat ${OUTNAME}.mt.1n.finalSnps.vcf | grep "PASS\|^#" > ${OUTNAME}.mt.1n.finalSnps_grepas.vcf

GenomeAnalysisTK.jar -T ValidateVariants \
-R $TARGETref \
--variant ${OUTNAME}.cp.1n.finalSnps.vcf

GenomeAnalysisTK.jar -T ValidateVariants \
-R $TARGETref \
--variant ${OUTNAME}.mt.1n.finalSnps.vcf
echo ''
echo 'Outfiles: '
echo '-> ${OUTNAME}.cp.finalSnps.vcf '
VCFGREPFILEcp=${OUTNAME}.cp.1n.finalSnps_grepas.vcf
VCFFILEcp=${OUTNAME}.cp.1n.finalSnps.vcf

echo '-> ${OUTNAME}.mt.finalSnps.vcf '
VCFGREPFILEmt=${OUTNAME}.mt.1n.finalSnps_grepas.vcf
VCFFILEmt=${OUTNAME}.mt.1n.finalSnps.vcf

echo ''
echo '###################################################################################'
echo "Get coverage"

bamToBed -i SWA1vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > SWA1vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed
bamToBed -i SWA2vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > SWA2vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed
bamToBed -i SWA3vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > SWA3vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed

# NC_013991.2 84 = cp
# NC_016740.1 40 = mt

grep "NC_016740" SWA1vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed > SWA1vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed
grep "NC_016740" SWA2vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed > SWA2vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed
grep "NC_016740" SWA3vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed > SWA3vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed

grep "NC_013991" SWA1vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed > SWA1vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed
grep "NC_013991" SWA2vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed > SWA2vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed
grep "NC_013991" SWA3vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed > SWA3vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed

# mt
samtools depth -b SWA1vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed SWA1vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > SWA1vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt
samtools depth -b SWA2vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed SWA2vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > SWA2vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt
samtools depth -b SWA3vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed SWA3vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > SWA3vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt

Rscript /scratch/nine/cov_musagenesAug14/coverage.R SWA1vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt SWA1vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup_cov.txt
Rscript /scratch/nine/cov_musagenesAug14/coverage.R SWA2vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt SWA2vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup_cov.txt
Rscript /scratch/nine/cov_musagenesAug14/coverage.R SWA3vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt SWA3vs124plgenes.mt.coordSorted_RGh.raln.recal3.sorted.Srmdup_cov.txt

# cp
samtools depth -b SWA1vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed SWA1vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > SWA1vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt
samtools depth -b SWA2vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed SWA2vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > SWA2vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt
samtools depth -b SWA3vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed SWA3vs124plgenes.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > SWA3vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt

Rscript /scratch/nine/cov_musagenesAug14/coverage.R SWA1vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt SWA1vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup_cov.txt
Rscript /scratch/nine/cov_musagenesAug14/coverage.R SWA2vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt SWA2vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup_cov.txt
Rscript /scratch/nine/cov_musagenesAug14/coverage.R SWA3vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt SWA3vs124plgenes.cp.coordSorted_RGh.raln.recal3.sorted.Srmdup_cov.txt

echo '########################## DONE ###################################################'