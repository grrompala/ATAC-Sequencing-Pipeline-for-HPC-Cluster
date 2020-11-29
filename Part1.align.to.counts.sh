# Notes on things needed to make this work
	# MutliQC code
	# installed multiqc locally with pip --user



# Necessary modules
ml fastqc trim_galore star samtools java picard python py_packages

# Working Directory
home=/sc/arion/projects/hmDNAmap/Dracheva.Heroin.Celltype.Project/ATACseq/Batch.A
output=/sc/arion/scratch/rompag01
cd ${output}
mkdir  STAR.output  MTreads processing.counts STAR.report duplicate.report fastqc trim.report
# STAR index folder
index=/sc/arion/projects/hmDNAmap/STAR.indices/STAR_hg19.2020
cd ${home}

for sample_folder in *GLU

do

# Grab Sample ID from file name
SAMPLE=$(echo ${sample_folder} | sed "s|Sample_pool-[A-Z]-||")
echo ${SAMPLE}

# Run FastQC	
cd ${output}/fastqc
fastqc -f fastq -o ${output}/fastqc ${home}/${sample_folder}/fastq/*.R1.fastq.gz
mv *${SAMPLE}*R1*html ${SAMPLE}.R1.html
mv *${SAMPLE}*R1*zip ${SAMPLE}.R1._fastqc.zip
fastqc -f fastq -o ${output}/fastqc ${home}/${sample_folder}/fastq/*.R2.fastq.gz
mv *${SAMPLE}*R2*html ${SAMPLE}.R2.html
mv *${SAMPLE}*R2*zip ${SAMPLE}.R2._fastqc.zip

# Run Trim Galore for Adapter trimming
# If space is an issue, can re-route trimmed files to a scratch folder or delete after using

cd ${home}

cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -q 20 --minimum-length 36 -j 8 \
-o ${output}/${SAMPLE}.R1.trimmed.fastq.gz -p ${output}/${SAMPLE}.R2.trimmed.fastq.gz \
${sample_folder}/fastq/*.R1.fastq.gz ${sample_folder}/fastq/*.R2.fastq.gz > ${output}/trim.report/${SAMPLE}.cutadapt.report.txt


# Genomic Alignment in Star
cd ${output}


STAR --genomeDir ${index} \
	 --runThreadN 8 \
	 --readFilesIn ${SAMPLE}.R1.trimmed.fastq.gz ${SAMPLE}.R2.trimmed.fastq.gz \
     --outFileNamePrefix ${output}/STAR.output/${SAMPLE}. \
	 --outSAMtype BAM SortedByCoordinate \
	 --readFilesCommand zcat

cd ${output}/STAR.output

# Move STAR output to folder
mv ${SAMPLE}*Log* ${SAMPLE}*tab ${output}/STAR.report
 
 # Grab mitochondria/chromosome read stats 
 samtools idxstats ${SAMPLE}.Aligned.sortedByCoord.out.bam > ${output}/MTreads/${SAMPLE}.txt
 
# Remove mitochondrial and multimapping reads from BAMs

 samtools view -h ${SAMPLE}.Aligned.sortedByCoord.out.bam | grep -v MT | samtools sort -O bam -o ${SAMPLE}_No.Mito.bam
 samtools view -bq 1  ${SAMPLE}_No.Mito.bam > ${SAMPLE}_No.Mito.Multi.bam
 # samtools view -b -F 30 ${SAMPLE}_No.Mito.bam > ${SAMPLE}.sorted.bam
 
 # Extracting read counts after each filtering step
 samtools view -c ${SAMPLE}.Aligned.sortedByCoord.out.bam > ${SAMPLE}.all.C
 samtools view -c ${SAMPLE}_No.Mito.bam > ${SAMPLE}.no.mito.C
 samtools view -c ${SAMPLE}_No.Mito.Multi.bam > ${SAMPLE}.no.mito.multi.C


# Remove duplicates

java -XX:ParallelGCThreads=8 \
     -jar ${PICARD} MarkDuplicates \
	 I=${SAMPLE}_No.Mito.Multi.bam \
	 O=${SAMPLE}.nodup.bam \
	 REMOVE_DUPLICATES=true \
	 METRICS_FILE=${output}/duplicate.report/${SAMPLE}.picard.metrics

samtools view -c ${SAMPLE}.nodup.bam > ${SAMPLE}.no.dup.C

# Shift Tn5
samtools index ${SAMPLE}.nodup.bam
alignmentSieve --numberOfProcessors 8 -b ${SAMPLE}.nodup.bam -o ${SAMPLE}.final.bam --ATACshift

cat ${SAMPLE}.all.C ${SAMPLE}.no.mito.C ${SAMPLE}.no.mito.multi.C ${SAMPLE}.no.dup.C > ${output}/processing.counts/${SAMPLE}.Processing.Counts

rm ${SAMPLE}.nodup.bam ${SAMPLE}.Aligned*bam *${SAMPLE}_No.Mito.bam ${SAMPLE}_No.Mito.Multi.bam
rm ${output}/STAR.output/${SAMPLE}*C
rm ${output}/STAR.output/*bai
rm ${output}/${SAMPLE}*trimmed*

done

cd ${BAM_dir}

multiqc STAR.report duplicate.report trim.report fastqc
