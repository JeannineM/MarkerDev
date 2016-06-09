#!/usr/bin/env bash
#usage time ./trim_align_gvcf.sh pairs_input.XX > trim_align_gvcf_XX.$(date +%d-%m).log.txt 2> trim_align_gvcf_XX.$(date +%d-%m).err.txt
# START in /data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/
#version March 2016
# 
#input files are the raw demultiplexed fastq.gz files
#say how many samples to be processed in pairs_input.XX


INPUT=$1 															# = pairs_input.XX
S=$(sed 's/.*\.//' <<< $INPUT) 										# extract the '.XX' from INPUT
OUTFOLD=trimmed_Mrz-16_${S}
TARGETref=refseq_gene236.fa

cd /data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/
mkdir $OUTFOLD

echo "###################################################################################"
echo '########################### trimmomatric ##########################################'
echo 'trimmed: T3PE2_myP_rc__anchored_MINL100_AVGQ30_R1R2'

cat $INPUT | while read line; \
do ID=$(sed 's/_.*//' <<< $line); echo ${ID}; \

trimmomatic-0.33.jar PE -phred33 \
$line \
${OUTFOLD}/${ID}.R1.trim_PE.fastq ${OUTFOLD}/${ID}.R1.trim_sing.fastq \
${OUTFOLD}/${ID}.R2.trim_PE.fastq ${OUTFOLD}/${ID}.R2.trim_sing.fastq \
ILLUMINACLIP:/data/scratch/btw680/analyses/reseq_ana/all_files/adapters/TruSeq3-PE-2.fa:2:30:10:1:TRUE ILLUMINACLIP:/data/scratch/btw680/analyses/reseq_ana/all_files/adapters/myBPrimer300.fa:2:30:10:1:TRUE \
ILLUMINACLIP:/data/scratch/btw680/analyses/reseq_ana/all_files/adapters/myBPrimer300_rc.fa:2:30:10:1:TRUE \
MINLEN:100 AVGQUAL:30 > ${OUTFOLD}/${ID}_trimmo_${S}.log.txt 2>&1; \

done


echo '##### trimming results #################'
cd ${OUTFOLD}
grep "Input Read Pairs" *_trimmo_${S}.log.txt > trimmo_${S}.result_$(date +%d-%m).tab    



cd /data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/${OUTFOLD}

cp /data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/all336/refseq_gene236.* ./
TARGETref=refseq_gene236.fa

cat ../$INPUT | while read line; \
do ID=$(sed 's/_.*//' <<< $line); \

LEFT=${ID}.R1.trim_PE.fastq; RIGHT=${ID}.R2.trim_PE.fastq; POPID=$(sed 's/-.*//' <<< $ID); OUTNAME=${ID}.vs236; \
echo ''; \
echo '###################################################################################'; \
echo '##### align PE reads to reference #################'; \
echo "@ started" $ID $POPID $OUTNAME $(date +%d-%m) $(date +%H:%M); \
alignReads.pl --left $LEFT --right $RIGHT --seqType fq --target $TARGETref --retain_SAM_files --output ${OUTNAME} --aligner bowtie2 -- -p 16 --dovetail --very-sensitive-local > ${ID}_realign_$(date +%d-%m).log.txt 2>&1; \ 

cd ${OUTNAME}; \
cp /data/scratch/btw680/analyses/reseq_ana/all_files/allRAW/all336/refseq_gene236.* ./ ; \
samtools view -hS -t target.fa.fai ${OUTNAME}.coordSorted.sam -o ${OUTNAME}.coordSorted_header.sam; \

AddOrReplaceReadGroups.jar I=${OUTNAME}.coordSorted_header.sam O=${OUTNAME}.coordSorted_RGh.sam SO=coordinate ID=${ID} PL=ILLUMINA SM=${ID} LB=L001 PU=${POPID} CREATE_INDEX=true; \

samtools view -bhS ${OUTNAME}.coordSorted_RGh.sam > ${OUTNAME}.coordSorted_RGh.bam; \
samtools index ${OUTNAME}.coordSorted_RGh.bam; \    

echo "###################################################################################"; \
echo "Realign reads around indels"; \

GenomeAnalysisTK.jar -T RealignerTargetCreator -R $TARGETref -dt NONE -I ${OUTNAME}.coordSorted_RGh.bam -o ${OUTNAME}.coordSorted_RGh.intervals; \
GenomeAnalysisTK.jar -T IndelRealigner -R $TARGETref -dt NONE -I ${OUTNAME}.coordSorted_RGh.bam -targetIntervals ${OUTNAME}.coordSorted_RGh.intervals -o ${OUTNAME}.coordSorted_RGh.realn.bam; \
samtools index ${OUTNAME}.coordSorted_RGh.realn.bam; \

samtools sort ${OUTNAME}.coordSorted_RGh.realn.bam ${OUTNAME}.coordSorted_RGh.realn.sorted; \
samtools index ${OUTNAME}.coordSorted_RGh.realn.sorted.bam; \

echo "###################################################################################"; \
echo "Recalibrate the base quality score (for multiple files)"; \
echo "First for BQSR of ${OUTNAME}.coordSorted_RGh.raln.bam"; \

GenomeAnalysisTK.jar -T UnifiedGenotyper -dt NONE -I ${OUTNAME}.coordSorted_RGh.realn.sorted.bam -R $TARGETref -glm BOTH -o ${OUTNAME}.coordSorted_RGh.realn.sorted.raw.vcf; \
GenomeAnalysisTK.jar -T BaseRecalibrator -dt NONE -I ${OUTNAME}.coordSorted_RGh.realn.sorted.bam -R $TARGETref -knownSites ${OUTNAME}.coordSorted_RGh.realn.sorted.raw.vcf --maximum_cycle_value 800 -o ${OUTNAME}.coordSorted_RGh.realn.sorted.recal_data.before.grp; \
GenomeAnalysisTK.jar -T PrintReads -R $TARGETref -dt NONE -I ${OUTNAME}.coordSorted_RGh.realn.sorted.bam -BQSR ${OUTNAME}.coordSorted_RGh.realn.sorted.recal_data.before.grp -o ${OUTNAME}.coordSorted_RGh.realn.sorted.recal.bam; \

echo ""; \
echo "Second round for BQSR of ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3.bam"; \

GenomeAnalysisTK.jar -T UnifiedGenotyper -R $TARGETref -dt NONE -I ${OUTNAME}.coordSorted_RGh.realn.sorted.recal.bam -glm BOTH -o ${OUTNAME}.coordSorted_RGh.realn.sorted.recal.raw.vcf; \
GenomeAnalysisTK.jar -T BaseRecalibrator -R $TARGETref -dt NONE -I ${OUTNAME}.coordSorted_RGh.realn.sorted.recal.bam -knownSites ${OUTNAME}.coordSorted_RGh.realn.sorted.recal.raw.vcf --maximum_cycle_value 800  -o ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2_data.before.grp; \
GenomeAnalysisTK.jar -T PrintReads -R $TARGETref -dt NONE -I ${OUTNAME}.coordSorted_RGh.realn.sorted.recal.bam -BQSR ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2_data.before.grp -o ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2.bam; \
GenomeAnalysisTK.jar -T BaseRecalibrator -R $TARGETref -dt NONE -I ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2.bam -knownSites ${OUTNAME}.coordSorted_RGh.realn.sorted.recal.raw.vcf --maximum_cycle_value 800 -o ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2_data.after.grp; \

echo ""; \
echo "Third round for BQSR of ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3.bam"; \

GenomeAnalysisTK.jar -T UnifiedGenotyper -R $TARGETref -dt NONE -I ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2.bam -glm BOTH -o ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2.raw.vcf; \
GenomeAnalysisTK.jar -T PrintReads -R $TARGETref -dt NONE -I ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2.bam -BQSR ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2_data.after.grp -o ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3.bam; \
GenomeAnalysisTK.jar -T BaseRecalibrator -R $TARGETref -dt NONE -I ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3.bam -knownSites ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2.raw.vcf --maximum_cycle_value 800 -o ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3_data.after.grp; \

echo "plot the BQSR comparison before and after"; \

GenomeAnalysisTK.jar -T AnalyzeCovariates \
-R $TARGETref \
-dt NONE \
-before ${OUTNAME}.coordSorted_RGh.realn.sorted.recal_data.before.grp \
-BQSR ${OUTNAME}.coordSorted_RGh.realn.sorted.recal2_data.before.grp \
-after ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3_data.after.grp \
-csv ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3_BQSR.csv \
-plots ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3_BQSR.pdf; \

echo "##################################################################################"; \
echo "Get coverage"; \

bamToBed -i ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3.bam > ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3.bed; \
samtools depth -b ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3.bed ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3.bam > ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3_samdepth.txt; \
Rscript /scratch/nine/cov_musagenesAug14/coverage.R ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3_samdepth.txt ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3_cov.txt; \

echo "##################################################################################"; \
echo "Get stats for read alignments"; \
samtools flagstat ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3.bam > out.bamstats.coord.txt ; \

echo "##################################################################################"; \
echo '# Gatk: call all whole set of single-sample gVCFs'; \

# .. -ploidy

GenomeAnalysisTK.jar -T HaplotypeCaller \
-R $TARGETref \
-I ${OUTNAME}.coordSorted_RGh.realn.sorted.recal3.bam \
-ERC GVCF \
-variant_index_type LINEAR -variant_index_parameter 128000 \
-bamout ${OUTNAME}.coordSorted_RGh.sorted.HC.bam \
-dt NONE \
-gt_mode DISCOVERY \
-stand_call_conf 30 \
-stand_emit_conf 20 \
-nda \
--max_alternate_alleles 22 \
-A Coverage -A MappingQualityZero -A FisherStrand -A QualByDepth -A DepthPerAlleleBySample \
-A ChromosomeCounts -A VariantType -A DepthPerSampleHC -A AlleleCountBySample -A StrandBiasBySample -A GenotypeSummaries \
-o ${ID}.gatkHC.raw.snps.indels_scc30_sec20.g.vcf & \

echo "##################################################################################"; \
echo "DONE" $ID; \
cd ../ ; \
done

echo "DONE:" $(date +%d-%m) $(date +%H:%M)