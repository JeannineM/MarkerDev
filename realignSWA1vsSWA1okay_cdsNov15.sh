#!/usr/bin/env bash
#usage time ./realign.sh > screen_$(date +%d-%m).log.txt 2>&1 &

### rename dict file for reference line 66

OUTDIR=/data/scratch/btw680/analyses/realign_WHtranscriptome
cd $OUTDIR

### refseqs from evigene -> main transcripts

TARGETdir=/data/scratch/btw680/analyses/trinBlast/inputreads/oksets_27June/SWA1_ren.okay_ren.cds
#TARGETdir=/data/scratch/btw680/analyses/trinBlast/inputreads/oksets_27June/SWA2_ren.okay_ren.cds
#TARGETdir=/data/scratch/btw680/analyses/trinBlast/inputreads/oksets_27June/SWA3_ren.okay_ren.cds
#TARGETdir=/data/scratch/btw680/analyses//evrogene/27-6-2014/okayset/SWA4_ren2.okay_ren.cds

#~/software/scripts/singleline.pl $TARGETdir > SWA1_ren.okay_ren.sl.cds

TARGETref=SWA1_ren.okay_ren.sl.cds
#TARGETref=SWA2_ren.okay_ren.sl.cds
#TARGETref=SWA3_ren.okay_ren.sl.cds
#TARGETref=SWA4_ren2.okay_ren.sl.cds

samtools faidx $TARGETref
#cp $TARGETref SWA1_ren.okay_ren.sl.fai

OUTNAME=SWA1vsSWA1okay_cds
LEFT=/data/home/btw680/archive/Jeannine_analyses/Illu_readsAprNew2014/SWA1_R1_trimcut6_stillpaired.fastq
RIGHT=/data/home/btw680/archive/Jeannine_analyses/Illu_readsAprNew2014/SWA1_R2_trimcut6_stillpaired.fastq

#OUTNAME=SWA2vsSWA2okay_cds
#LEFT=/data/home/btw680/archive/Jeannine_analyses/Illu_readsAprNew2014/SWA2_R1_trimcut11_stillpaired.fastq 
#RIGHT=/data/home/btw680/archive/Jeannine_analyses/Illu_readsAprNew2014/SWA2_R2_trimcut12_stillpaired.fastq

#OUTNAME=SWA3vsSWA3okay_cds
#LEFT=/data/home/btw680/archive/Jeannine_analyses/Illu_readsAprNew2014/SWA3_R1_trimcut12_stillpaired.fastq
#RIGHT=/data/home/btw680/archive/Jeannine_analyses/Illu_readsAprNew2014/SWA3_R2_trimcut12_stillpaired.fastq

#OUTNAME=SWA4vsSWA4okay_cds
#SINGLE=/data/home/btw680/archive/Jeannine_analyses/Illu_readsAprNew2014/SWAall_singlets.fastq
#SINGLE=/data/home/btw680/archive/Jeannine_analyses/Illu_readsAprNew2014/454_trimmedMay2014.fastq

#cd $OUTDIR
#mkdir $OUTNAME
#cd $OUTDIR/$OUTNAME

# echo "###################################################################################"
# echo "SWA1 read align to $TARGETref"
# cd $OUTDIR
# #~/software/scripts/singleline.pl SWA1_ren.okay_ren.cds > SWA1_ren.okay_ren.sl.cds
# #samtools faidx $TARGETref
# 
# time alignReads.pl --left $LEFT --right $RIGHT --seqType fq --target $TARGETref --retain_SAM_files --output $OUTNAME --aligner bowtie2 -- -p 16 --very-sensitive > ${OUTNAME}_$(date +%d-%m).log.txt 2>&1
# 
# cd $OUTDIR/$OUTNAME
# rm *.finished
# 
# echo "###################################################################################"
# echo "Output stats on ${OUTNAME}.coordSorted.sam"
# cd $OUTDIR/$OUTNAME
# 
# SAM_nameSorted_to_uniq_count_stats.pl ${OUTNAME}.nameSorted.sam > out.SAMname.stats.txt
# samstats.pl ${OUTNAME}.nameSorted.sam > out.samstats.coord.txt
# 
# echo "###################################################################################"
# echo "Prepare reference for Gatk input"
# #####-> 1. the reference file
cd $OUTDIR/$OUTNAME
cp ../$TARGETref ./
# 
mv $TARGETref ${TARGETref}.fa
CreateSequenceDictionary.jar R=${TARGETref}.fa O=${TARGETref}.dict
cp ${TARGETref}.dict SWA1_ren.okay_ren.sl.dict
# 
# samtools faidx ${TARGETref}.fa
# 
# echo "###################################################################################"
# echo "Prepare SWA1 sam/bam for Gatk input"
# #->the sam / bam input file ## from the trinity alignReads.pl output:
 cd $OUTDIR/$OUTNAME
# 
samtools view -hS -t ${TARGETref}.fa.fai ${OUTNAME}.coordSorted.sam -o ${OUTNAME}.coordSorted_header.sam

AddOrReplaceReadGroups.jar I=${OUTNAME}.coordSorted_header.sam O=${OUTNAME}.coordSorted_RGh.sam SO=coordinate ID=SWA1 PL=ILLUMINA SM=SWA1 LB=LIB-SWA1-1 PI=202 PU=SWA1 CREATE_INDEX=true
#AddOrReplaceReadGroups.jar I=${OUTNAME}.coordSorted_header.sam O=${OUTNAME}.coordSorted_RGh.sam SO=coordinate ID=SWA2 PL=ILLUMINA SM=SWA2 LB=LIB-SWA2-1 PI=240 PU=SWA2 CREATE_INDEX=true
#AddOrReplaceReadGroups.jar I=${OUTNAME}.coordSorted_header.sam O=${OUTNAME}.coordSorted_RGh.sam SO=coordinate ID=SWA3 PL=ILLUMINA SM=SWA3 LB=LIB-SWA3-1 PI=228 PU=SWA3 CREATE_INDEX=true


ValidateSamFile.jar INPUT=${OUTNAME}.coordSorted_RGh.sam 

samtools view -bhS ${OUTNAME}.coordSorted_RGh.sam > ${OUTNAME}.coordSorted_RGh.bam
samtools index ${OUTNAME}.coordSorted_RGh.bam

# 
echo "###################################################################################"
echo "Realign reads around indels"
# 
# GenomeAnalysisTK.jar -T RealignerTargetCreator -R ${TARGETref}.fa -I ${OUTNAME}.coordSorted_RGh.bam -o ${OUTNAME}.coordSorted_RGh.intervals
# GenomeAnalysisTK.jar -T IndelRealigner -R ${TARGETref}.fa -I ${OUTNAME}.coordSorted_RGh.bam -targetIntervals ${OUTNAME}.coordSorted_RGh.intervals -o ${OUTNAME}.coordSorted_RGh.raln.bam
samtools index ${OUTNAME}.coordSorted_RGh.raln.bam 

echo "###################################################################################"
echo "Recalibrate the base quality score (for multiple files)"
echo "first for BQSR of ${OUTNAME}.coordSorted_RGh.raln.bam"

GenomeAnalysisTK.jar -T UnifiedGenotyper -I ${OUTNAME}.coordSorted_RGh.raln.bam -R ${TARGETref}.fa -glm BOTH -o ${OUTNAME}.coordSorted_RGh.raln.raw.vcf
GenomeAnalysisTK.jar -T BaseRecalibrator -I ${OUTNAME}.coordSorted_RGh.raln.bam -R ${TARGETref}.fa -knownSites ${OUTNAME}.coordSorted_RGh.raln.raw.vcf --maximum_cycle_value 800 -o ${OUTNAME}.coordSorted_RGh.raln.recal_data.before.grp
GenomeAnalysisTK.jar -T PrintReads -R ${TARGETref}.fa -I ${OUTNAME}.coordSorted_RGh.raln.bam -BQSR ${OUTNAME}.coordSorted_RGh.raln.recal_data.before.grp -o ${OUTNAME}.coordSorted_RGh.raln.recal.bam

echo ""
echo "second round for BQSR of ${OUTNAME}.coordSorted_RGh.raln.bam"

GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${TARGETref}.fa -I ${OUTNAME}.coordSorted_RGh.raln.recal.bam -glm BOTH -o ${OUTNAME}.coordSorted_RGh.raln.recal.raw.vcf
GenomeAnalysisTK.jar -T BaseRecalibrator -R ${TARGETref}.fa -I ${OUTNAME}.coordSorted_RGh.raln.recal.bam -knownSites ${OUTNAME}.coordSorted_RGh.raln.recal.raw.vcf --maximum_cycle_value 800  -o ${OUTNAME}.coordSorted_RGh.raln.recal2_data.before.grp
GenomeAnalysisTK.jar -T PrintReads -R ${TARGETref}.fa -I ${OUTNAME}.coordSorted_RGh.raln.recal.bam -BQSR ${OUTNAME}.coordSorted_RGh.raln.recal2_data.before.grp -o ${OUTNAME}.coordSorted_RGh.raln.recal2.bam
GenomeAnalysisTK.jar -T BaseRecalibrator -R ${TARGETref}.fa -I ${OUTNAME}.coordSorted_RGh.raln.recal2.bam -knownSites ${OUTNAME}.coordSorted_RGh.raln.recal.raw.vcf --maximum_cycle_value 800 -o ${OUTNAME}.coordSorted_RGh.raln.recal2_data.after.grp

echo ""
echo "third round for BQSR of ${OUTNAME}.coordSorted_RGh.raln.bam"

GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${TARGETref}.fa -I ${OUTNAME}.coordSorted_RGh.raln.recal2.bam -glm BOTH -o ${OUTNAME}.coordSorted_RGh.raln.recal2.raw.vcf
GenomeAnalysisTK.jar -T PrintReads -R ${TARGETref}.fa -I ${OUTNAME}.coordSorted_RGh.raln.recal2.bam -BQSR ${OUTNAME}.coordSorted_RGh.raln.recal2_data.after.grp -o ${OUTNAME}.coordSorted_RGh.raln.recal3.bam
GenomeAnalysisTK.jar -T BaseRecalibrator -R ${TARGETref}.fa -I ${OUTNAME}.coordSorted_RGh.raln.recal3.bam -knownSites ${OUTNAME}.coordSorted_RGh.raln.recal2.raw.vcf --maximum_cycle_value 800 -o ${OUTNAME}.coordSorted_RGh.raln.recal3_data.after.grp

echo "plot the BQSR comparison before and after"
cd $OUTDIR/$OUTNAME

GenomeAnalysisTK.jar -T AnalyzeCovariates \
-R ${TARGETref}.fa \
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