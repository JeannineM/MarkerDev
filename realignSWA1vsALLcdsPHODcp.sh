#!/usr/bin/env bash
#usage bash realign.sh > screen.log 2>&1
# 

TARGETref=PhoD_cp_nt_cds_ren.fa 
LEFT=/data/home/btw680/archive/Jeannine_analyses/Illu_readsAprNew2014/SWA1_R1_trimcut6_stillpaired.fastq
RIGHT=/data/home/btw680/archive/Jeannine_analyses/Illu_readsAprNew2014/SWA1_R2_trimcut6_stillpaired.fastq
OUTNAME=SWA1vsALLcdsPHODcp

cd /scratch/nine/cov_musagenesAug14/cpalign/
mkdir SWAvsALLcdsPHODcp
cd /scratch/nine/cov_musagenesAug14/cpalign/SWAvsALLcdsPHODcp
cp /scratch/nine/cov_musagenesAug14/cpalign/$TARGETref ./

mkdir $OUTNAME
cd /scratch/nine/cov_musagenesAug14/cpalign/SWAvsALLcdsPHODcp/$OUTNAME


echo "###################################################################################"
echo "SWA1read align to $TARGETref"
cd /scratch/nine/cov_musagenesAug14/cpalign/SWAvsALLcdsPHODcp

time alignReads.pl --left $LEFT --right $RIGHT --seqType fq --target $TARGETref --retain_SAM_files --output $OUTNAME --aligner bowtie2 -- -p 16 --very-sensitive > ${OUTNAME}_$(date +%d-%m).log.txt 2>&1

cd /scratch/nine/cov_musagenesAug14/cpalign/SWAvsALLcdsPHODcp/$OUTNAME
rm *.finished

echo "###################################################################################"
echo "Output stats on ${OUTNAME}.coordSorted.sam"
cd /scratch/nine/cov_musagenesAug14/cpalign/SWAvsALLcdsPHODcp/$OUTNAME

SAM_nameSorted_to_uniq_count_stats.pl ${OUTNAME}.nameSorted.sam > out.SAMname.stats.txt
samstats.pl ${OUTNAME}.nameSorted.sam > out.samstats.coord.txt

echo "###################################################################################"
echo "Prepare reference for Gatk input"
#####-> 1. the reference file
cd /scratch/nine/cov_musagenesAug14/cpalign/SWAvsALLcdsPHODcp
#cp $TARGETref ./

CreateSequenceDictionary.jar R=$TARGETref O=${TARGETref}.dict
cp ${TARGETref}.dict  PhoD_cp_nt_cds_ren.dict

samtools faidx $TARGETref
# 
echo "###################################################################################"
echo "Prepare SWA1 sam/bam for Gatk input"
#->the sam / bam input file ## from the trinity alignReads.pl output:
cd /scratch/nine/cov_musagenesAug14/cpalign/SWAvsALLcdsPHODcp/$OUTNAME

samtools view -hS -t ../${TARGETref}.fai ${OUTNAME}.coordSorted.sam -o ${OUTNAME}.coordSorted_header.sam

AddOrReplaceReadGroups.jar I=${OUTNAME}.coordSorted_header.sam O=${OUTNAME}.coordSorted_RGh.sam SO=coordinate ID=FC1.L1 PL=ILLUMINA SM=SWA1 LB=LIB-SWA1-1 PI=250 PU=SWA1 CREATE_INDEX=true

ValidateSamFile.jar INPUT=${OUTNAME}.coordSorted_RGh.sam 

samtools view -bhS ${OUTNAME}.coordSorted_RGh.sam > ${OUTNAME}.coordSorted_RGh.bam
samtools index ${OUTNAME}.coordSorted_RGh.bam

echo "###################################################################################"
echo "Realign reads around indels"

GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../$TARGETref -I ${OUTNAME}.coordSorted_RGh.bam -o ${OUTNAME}.coordSorted_RGh.intervals
GenomeAnalysisTK.jar -T IndelRealigner -R ../$TARGETref -I ${OUTNAME}.coordSorted_RGh.bam -targetIntervals ${OUTNAME}.coordSorted_RGh.intervals -o ${OUTNAME}.coordSorted_RGh.raln.bam
samtools index ${OUTNAME}.coordSorted_RGh.raln.bam 

echo "###################################################################################"
echo "Recalibrate the base quality score (for multiple files)"
echo "first for BQSR of ${OUTNAME}.coordSorted_RGh.raln.bam"

GenomeAnalysisTK.jar -T UnifiedGenotyper -I ${OUTNAME}.coordSorted_RGh.raln.bam -R ../$TARGETref -glm BOTH -o ${OUTNAME}.coordSorted_RGh.raln.raw.vcf
GenomeAnalysisTK.jar -T BaseRecalibrator -I ${OUTNAME}.coordSorted_RGh.raln.bam -R ../$TARGETref -knownSites ${OUTNAME}.coordSorted_RGh.raln.raw.vcf -o ${OUTNAME}.coordSorted_RGh.raln.recal_data.before.grp
GenomeAnalysisTK.jar -T PrintReads -R ../$TARGETref -I ${OUTNAME}.coordSorted_RGh.raln.bam -BQSR ${OUTNAME}.coordSorted_RGh.raln.recal_data.before.grp -o ${OUTNAME}.coordSorted_RGh.raln.recal.bam

echo ""
echo "second round for BQSR of ${OUTNAME}.coordSorted_RGh.raln.bam"

GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../$TARGETref -I ${OUTNAME}.coordSorted_RGh.raln.recal.bam -glm BOTH -o ${OUTNAME}.coordSorted_RGh.raln.recal.raw.vcf
GenomeAnalysisTK.jar -T BaseRecalibrator -R ../$TARGETref -I ${OUTNAME}.coordSorted_RGh.raln.recal.bam -knownSites ${OUTNAME}.coordSorted_RGh.raln.recal.raw.vcf -o ${OUTNAME}.coordSorted_RGh.raln.recal2_data.before.grp
GenomeAnalysisTK.jar -T PrintReads -R ../$TARGETref -I ${OUTNAME}.coordSorted_RGh.raln.recal.bam -BQSR ${OUTNAME}.coordSorted_RGh.raln.recal2_data.before.grp -o ${OUTNAME}.coordSorted_RGh.raln.recal2.bam
GenomeAnalysisTK.jar -T BaseRecalibrator -R ../$TARGETref -I ${OUTNAME}.coordSorted_RGh.raln.recal2.bam -knownSites ${OUTNAME}.coordSorted_RGh.raln.recal.raw.vcf -o ${OUTNAME}.coordSorted_RGh.raln.recal2_data.after.grp

echo ""
echo "third round for BQSR of ${OUTNAME}.coordSorted_RGh.raln.bam"

GenomeAnalysisTK.jar -T UnifiedGenotyper -R ../$TARGETref -I ${OUTNAME}.coordSorted_RGh.raln.recal2.bam -glm BOTH -o ${OUTNAME}.coordSorted_RGh.raln.recal2.raw.vcf
GenomeAnalysisTK.jar -T PrintReads -R ../$TARGETref -I ${OUTNAME}.coordSorted_RGh.raln.recal2.bam -BQSR ${OUTNAME}.coordSorted_RGh.raln.recal2_data.after.grp -o ${OUTNAME}.coordSorted_RGh.raln.recal3.bam
GenomeAnalysisTK.jar -T BaseRecalibrator -R ../$TARGETref -I ${OUTNAME}.coordSorted_RGh.raln.recal3.bam -knownSites ${OUTNAME}.coordSorted_RGh.raln.recal2.raw.vcf -o ${OUTNAME}.coordSorted_RGh.raln.recal3_data.after.grp

echo "plot the BQSR comparison before and after"
cd /scratch/nine/cov_musagenesAug14/cpalign/SWAvsALLcdsPHODcp/$OUTNAME

GenomeAnalysisTK.jar -T AnalyzeCovariates \
-R ../$TARGETref \
-before ${OUTNAME}.coordSorted_RGh.raln.recal_data.before.grp \
-BQSR ${OUTNAME}.coordSorted_RGh.raln.recal2_data.before.grp \
-after ${OUTNAME}.coordSorted_RGh.raln.recal3_data.after.grp \
-csv ${OUTNAME}.coordSorted_RGh.raln.recal3_BQSR.csv \
-plots ${OUTNAME}.coordSorted_RGh.raln.recal3_BQSR.pdf

echo "##################################################################################"
echo "Coord sort & removal of non-unique aligned reads (rmdup only requires bam file)"
samtools sort ${OUTNAME}.coordSorted_RGh.raln.recal3.bam ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted
samtools index ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.bam
samtools rmdup -S ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.bam ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam
samtools index ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam 

#-> or masking of duplicate reads: it masks ~ 17% 
#MarkDuplicates.jar I=${OUTNAME}.coordSorted_RGh.raln.bam O=${OUTNAME}.coordSorted_RGh.raln.dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

echo "##################################################################################"
echo "Get coverage"

bamToBed -i ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed
samtools depth -b ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bed ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam > ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt
Rscript /scratch/nine/cov_musagenesAug14/coverage.R ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup_samdepth.txt ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup_cov.txt

echo "##################################################################################"
echo "Get stats"
echo "before:"
echo "> ${OUTNAME}.coordSorted.bam:"
samtools flagstat ${OUTNAME}.coordSorted.bam

echo "after:"
echo "> ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam:"
samtools flagstat ${OUTNAME}.coordSorted_RGh.raln.recal3.sorted.Srmdup.bam

echo "##################################################################################"
echo "DONE"