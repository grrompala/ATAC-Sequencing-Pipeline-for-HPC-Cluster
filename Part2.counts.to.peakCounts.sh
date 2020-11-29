## Starting with directory of BAM files
BAM_dir=/sc/arion/scratch/rompag01/STAR.output

mkdir genrich_peaks junk

cd ${BAM_dir}

#Genrich peak-calling
ml genrich
ml bedtools
ml gnuParallel
ml samtools

parallel 'samtools sort -@3 {} -o {.} -n' ::: *.bam
for file in *Sieve
do
mv ${file} ${file/.Sieve/.bam}
done
mkdir move
mv *Sieve.bam move
parallel 'Genrich -t {} -o ../genrich_peaks/{.}.genrich -j -D' ::: *.bam

# Merge all BAM files 
cd ${BAM_dir}

samtools merge -@8 MERGE.bam ${BAM_dir}/*bam 

Genrich -t MERGE.bam -o ../genrich_peaks/MERGE.genrich -v -j -D

bedtools intersect -v -a ../genrich_peaks/MERGE.genrich -b ../hg19-blacklist.bed > ../MERGED.genrich

mv ../genrich_peaks/MERGE.genich junk

# Run R script for consensus peak-calling
ml R

Rscript ../consensus.builder.R

echo "Starting FeatureCounts"

rm MERGE.bam

ml subread
featureCounts -F SAF -a ../Consensus_Peak.SAF -o ../Peak.counts.txt -p *.bam -Q 20


