module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 
java -XX:ParallelGCThreads=1 -Xmx2G -jar ${PICARD_HOME}/ReorderSam.jar \
  VALIDATION_STRINGENCY=SILENT \
  INPUT=alignment/brain/tophat/brain.sorted.bam \
  OUTPUT=alignment/brain/brain.sorted.bam \
  REFERENCE=/software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa

java -XX:ParallelGCThreads=1 -Xmx2G -jar ${PICARD_HOME}/ReorderSam.jar \
  VALIDATION_STRINGENCY=SILENT \
  INPUT=alignment/adrenal/tophat/adrenal.sorted.bam \
  OUTPUT=alignment/adrenal/adrenal.sorted.bam \
  REFERENCE=/software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa

#-------------------------------------------------------------------------------
# JOB: picard_sort_sam.brain
#-------------------------------------------------------------------------------
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/brain/brain.sorted.bam \
  OUTPUT=alignment/brain/brain.QueryNameSorted.bam \
  SORT_ORDER=queryname \
  MAX_RECORDS_IN_RAM=6750000


#-------------------------------------------------------------------------------
# JOB: picard_sort_sam.adrenal
#-------------------------------------------------------------------------------
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $PICARD_HOME/SortSam.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/adrenal/adrenal.sorted.bam \
  OUTPUT=alignment/adrenal/adrenal.QueryNameSorted.bam \
  SORT_ORDER=queryname \
  MAX_RECORDS_IN_RAM=6750000


#-------------------------------------------------------------------------------
# STEP: picard_mark_duplicates
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates.brain
#-------------------------------------------------------------------------------
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/brain/brain.sorted.bam \
  OUTPUT=alignment/brain/brain.sorted.mdup.bam \
  METRICS_FILE=alignment/brain/brain.sorted.mdup.metrics \
  MAX_RECORDS_IN_RAM=350000


#-------------------------------------------------------------------------------
# JOB: picard_mark_duplicates.adrenal
#-------------------------------------------------------------------------------
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $PICARD_HOME/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/adrenal/adrenal.sorted.bam \
  OUTPUT=alignment/adrenal/adrenal.sorted.mdup.bam \
  METRICS_FILE=alignment/adrenal/adrenal.sorted.mdup.metrics \
  MAX_RECORDS_IN_RAM=350000


#-------------------------------------------------------------------------------
# STEP: rnaseqc
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: rnaseqc
#-------------------------------------------------------------------------------
module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/bwa/0.7.10 mugqic/rnaseqc/1.1.7 && \
mkdir -p metrics/rnaseqRep && \
echo "Sample	BamFile	Note
brain	alignment/brain/brain.sorted.mdup.bam	RNAseq
adrenal	alignment/adrenal/adrenal.sorted.mdup.bam	RNAseq" \
  > alignment/rnaseqc.samples.txt && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $RNASEQC_JAR \
  -BWArRNA /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/ncrna_bwa_index/Homo_sapiens.GRCh37.Ensembl75.ncrna.fa \
  -n 1000 \
  -o metrics/rnaseqRep \
  -r /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  -s alignment/rnaseqc.samples.txt \
  -t /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Ensembl75.transcript_id.gtf \
  -ttype 2 && \
zip -r metrics/rnaseqRep.zip metrics/rnaseqRep


#-------------------------------------------------------------------------------
# STEP: wiggle
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: wiggle.brain.forward_strandspec
#-------------------------------------------------------------------------------
module load mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
samtools view -bh -F 256 -f 81 \
  alignment/brain/brain.sorted.mdup.bam \
  > alignment/brain/brain.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/brain/brain.sorted.mdup.bam \
  > alignment/brain/brain.sorted.mdup.tmp2.forward.bam && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  INPUT=alignment/brain/brain.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/brain/brain.sorted.mdup.tmp2.forward.bam \
  OUTPUT=alignment/brain/brain.sorted.mdup.forward.bam \
  MAX_RECORDS_IN_RAM=6750000 && \
rm alignment/brain/brain.sorted.mdup.tmp1.forward.bam alignment/brain/brain.sorted.mdup.tmp2.forward.bam


#-------------------------------------------------------------------------------
# JOB: wiggle.brain.reverse_strandspec
#-------------------------------------------------------------------------------
module load mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
mkdir -p tracks/brain tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/brain/brain.sorted.mdup.bam \
  > alignment/brain/brain.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/brain/brain.sorted.mdup.bam \
  > alignment/brain/brain.sorted.mdup.tmp2.reverse.bam && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  INPUT=alignment/brain/brain.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/brain/brain.sorted.mdup.tmp2.reverse.bam \
  OUTPUT=alignment/brain/brain.sorted.mdup.reverse.bam \
  MAX_RECORDS_IN_RAM=6750000 && \
rm alignment/brain/brain.sorted.mdup.tmp1.reverse.bam alignment/brain/brain.sorted.mdup.tmp2.reverse.bam


#-------------------------------------------------------------------------------
# JOB: wiggle.brain.forward
#-------------------------------------------------------------------------------
module load mugqic/samtools/0.1.19 mugqic/bedtools/2.17.0 mugqic/ucsc/20140212 && \
mkdir -p tracks/brain tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81 alignment/brain/brain.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/brain/brain.sorted.mdup.bam \
  -g /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa.fai \
  > tracks/brain/brain.forward.bedGraph && \
bedGraphToBigWig \
  tracks/brain/brain.forward.bedGraph \
  /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa.fai \
  tracks/bigWig/brain.forward.bw


#-------------------------------------------------------------------------------
# JOB: wiggle.brain.reverse
#-------------------------------------------------------------------------------
module load mugqic/samtools/0.1.19 mugqic/bedtools/2.17.0 mugqic/ucsc/20140212 && \
mkdir -p tracks/brain tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81 alignment/brain/brain.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/brain/brain.sorted.mdup.bam \
  -g /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa.fai \
  > tracks/brain/brain.reverse.bedGraph && \
bedGraphToBigWig \
  tracks/brain/brain.reverse.bedGraph \
  /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa.fai \
  tracks/bigWig/brain.reverse.bw


#-------------------------------------------------------------------------------
# JOB: wiggle.adrenal.forward_strandspec
#-------------------------------------------------------------------------------
module load mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
samtools view -bh -F 256 -f 81 \
  alignment/adrenal/adrenal.sorted.mdup.bam \
  > alignment/adrenal/adrenal.sorted.mdup.tmp1.forward.bam && \
samtools view -bh -F 256 -f 161 \
  alignment/adrenal/adrenal.sorted.mdup.bam \
  > alignment/adrenal/adrenal.sorted.mdup.tmp2.forward.bam && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  INPUT=alignment/adrenal/adrenal.sorted.mdup.tmp1.forward.bam \
  INPUT=alignment/adrenal/adrenal.sorted.mdup.tmp2.forward.bam \
  OUTPUT=alignment/adrenal/adrenal.sorted.mdup.forward.bam \
  MAX_RECORDS_IN_RAM=6750000 && \
rm alignment/adrenal/adrenal.sorted.mdup.tmp1.forward.bam alignment/adrenal/adrenal.sorted.mdup.tmp2.forward.bam


#-------------------------------------------------------------------------------
# JOB: wiggle.adrenal.reverse_strandspec
#-------------------------------------------------------------------------------
module load mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
mkdir -p tracks/adrenal tracks/bigWig && \
samtools view -bh -F 256 -f 97 \
  alignment/adrenal/adrenal.sorted.mdup.bam \
  > alignment/adrenal/adrenal.sorted.mdup.tmp1.reverse.bam && \
samtools view -bh -F 256 -f 145 \
  alignment/adrenal/adrenal.sorted.mdup.bam \
  > alignment/adrenal/adrenal.sorted.mdup.tmp2.reverse.bam && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $PICARD_HOME/MergeSamFiles.jar \
  VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true CREATE_INDEX=true \
  INPUT=alignment/adrenal/adrenal.sorted.mdup.tmp1.reverse.bam \
  INPUT=alignment/adrenal/adrenal.sorted.mdup.tmp2.reverse.bam \
  OUTPUT=alignment/adrenal/adrenal.sorted.mdup.reverse.bam \
  MAX_RECORDS_IN_RAM=6750000 && \
rm alignment/adrenal/adrenal.sorted.mdup.tmp1.reverse.bam alignment/adrenal/adrenal.sorted.mdup.tmp2.reverse.bam


#-------------------------------------------------------------------------------
# JOB: wiggle.adrenal.forward
#-------------------------------------------------------------------------------
module load mugqic/samtools/0.1.19 mugqic/bedtools/2.17.0 mugqic/ucsc/20140212 && \
mkdir -p tracks/adrenal tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81 alignment/adrenal/adrenal.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/adrenal/adrenal.sorted.mdup.bam \
  -g /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa.fai \
  > tracks/adrenal/adrenal.forward.bedGraph && \
bedGraphToBigWig \
  tracks/adrenal/adrenal.forward.bedGraph \
  /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa.fai \
  tracks/bigWig/adrenal.forward.bw


#-------------------------------------------------------------------------------
# JOB: wiggle.adrenal.reverse
#-------------------------------------------------------------------------------
module load mugqic/samtools/0.1.19 mugqic/bedtools/2.17.0 mugqic/ucsc/20140212 && \
mkdir -p tracks/adrenal tracks/bigWig && \
nmblines=$(samtools view -F 256 -f 81 alignment/adrenal/adrenal.sorted.mdup.bam | wc -l) && \
scalefactor=0$(echo "scale=2; 1 / ($nmblines / 10000000);" | bc) && \
genomeCoverageBed -bg -split -scale $scalefactor \
  -ibam alignment/adrenal/adrenal.sorted.mdup.bam \
  -g /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa.fai \
  > tracks/adrenal/adrenal.reverse.bedGraph && \
bedGraphToBigWig \
  tracks/adrenal/adrenal.reverse.bedGraph \
  /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa.fai \
  tracks/bigWig/adrenal.reverse.bw


#-------------------------------------------------------------------------------
# STEP: raw_counts
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: htseq_count.brain
#-------------------------------------------------------------------------------
echo "module load mugqic/samtools/0.1.19 mugqic/python/2.7.8 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/brain/brain.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Ensembl75.gtf \
  > raw_counts/brain.readcounts.csv" | qsub -V -m ae -M $JOB_MAIL -d . -j oe -N brain.htseq -l walltime=24:00:00 -q metaq -l nodes=1:ppn=2


#-------------------------------------------------------------------------------
# JOB: htseq_count.adrenal
#-------------------------------------------------------------------------------
echo "module load mugqic/samtools/0.1.19 mugqic/python/2.7.8 && \
mkdir -p raw_counts && \
samtools view -F 4 \
  alignment/adrenal/adrenal.QueryNameSorted.bam | \
htseq-count - \
  -m intersection-nonempty \
  --stranded=reverse \
  --format=sam \
  /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Ensembl75.gtf \
  > raw_counts/adrenal.readcounts.csv"  | qsub -V -m ae -M $JOB_MAIL -d . -j oe -N adrenal.htseq -l walltime=24:00:00 -q metaq -l nodes=1:ppn=2


#-------------------------------------------------------------------------------
# STEP: raw_counts_metrics
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: metrics.matrix
#-------------------------------------------------------------------------------
module load mugqic/mugqic_tools/2.0.2 && \
mkdir -p DGE && \
gtf2tmpMatrix.awk \
  /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Ensembl75.gtf \
  DGE/tmpMatrix.txt && \
HEAD='Gene\tSymbol' && \
for read_count_file in \
  raw_counts/brain.readcounts.csv \
  raw_counts/adrenal.readcounts.csv
do
  sort -k1,1 $read_count_file > DGE/tmpSort.txt && \
  join -1 1 -2 1 <(sort -k1,1 DGE/tmpMatrix.txt) DGE/tmpSort.txt > DGE/tmpMatrix.2.txt && \
  mv DGE/tmpMatrix.2.txt DGE/tmpMatrix.txt && \
  na=$(basename $read_count_file | cut -d. -f1) && \
  HEAD="$HEAD\t$na"
done && \
echo -e $HEAD | cat - DGE/tmpMatrix.txt | tr ' ' '\t' > DGE/rawCountMatrix.csv && \
rm DGE/tmpSort.txt DGE/tmpMatrix.txt


#-------------------------------------------------------------------------------
# JOB: metrics.wigzip
#-------------------------------------------------------------------------------
zip -r tracks.zip tracks/bigWig


#-------------------------------------------------------------------------------
# JOB: rpkm_saturation
#-------------------------------------------------------------------------------
module load mugqic/R_Bioconductor/3.1.2_3.0 mugqic/mugqic_tools/2.0.2 && \
mkdir -p metrics/saturation && \
Rscript $R_TOOLS/rpkmSaturation.R \
  DGE/rawCountMatrix.csv \
  /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Ensembl75.genes.length.tsv \
  raw_counts \
  metrics/saturation \
  2 \
  1 && \
zip -r metrics/saturation.zip metrics/saturation


#-------------------------------------------------------------------------------
# STEP: cufflinks
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: cufflinks.brain
#-------------------------------------------------------------------------------
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/brain && \
cufflinks -q  \
  --GTF-guide /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Ensembl75.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/brain \
  --num-threads 1 \
  alignment/brain/brain.sorted.mdup.bam


#-------------------------------------------------------------------------------
# JOB: cufflinks.adrenal
#-------------------------------------------------------------------------------
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/adrenal && \
cufflinks -q  \
  --GTF-guide /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Ensembl75.gtf \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/adrenal \
  --num-threads 1 \
  alignment/adrenal/adrenal.sorted.mdup.bam


#-------------------------------------------------------------------------------
# STEP: cuffmerge
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: cuffmerge
#-------------------------------------------------------------------------------
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/AllSamples && \
`cat > cufflinks/cuffmerge.samples.txt << END
cufflinks/brain/transcripts.gtf
cufflinks/adrenal/transcripts.gtf
END
  
` && \
cuffmerge  \
  --ref-gtf /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Ensembl75.gtf \
  --ref-sequence /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  -o cufflinks/AllSamples \
  --num-threads 1 \
  cufflinks/cuffmerge.samples.txt


#-------------------------------------------------------------------------------
# STEP: cuffquant
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: cuffquant.brain
#-------------------------------------------------------------------------------
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/brain && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/brain \
  --num-threads 1 \
  cufflinks/AllSamples/merged.gtf \
  alignment/brain/brain.sorted.mdup.bam


#-------------------------------------------------------------------------------
# JOB: cuffquant.adrenal
#-------------------------------------------------------------------------------
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cufflinks/adrenal && \
cuffquant -q  \
  --max-bundle-frags 1000000 \
  --library-type fr-firststrand \
  --output-dir cufflinks/adrenal \
  --num-threads 1 \
  cufflinks/AllSamples/merged.gtf \
  alignment/adrenal/adrenal.sorted.mdup.bam


#-------------------------------------------------------------------------------
# STEP: cuffdiff
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: cuffdiff.Contrast1
#-------------------------------------------------------------------------------
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffdiff/Contrast1 && \
cuffdiff -u \
  --frag-bias-correct /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/genome/Homo_sapiens.GRCh37.fa \
  --library-type fr-firststrand \
  --output-dir cuffdiff/Contrast1 \
  --num-threads 1 \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/brain/abundances.cxb \
  cufflinks/adrenal/abundances.cxb


#-------------------------------------------------------------------------------
# STEP: cuffnorm
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: cuffnorm
#-------------------------------------------------------------------------------
module load mugqic/cufflinks/2.2.1 && \
mkdir -p cuffnorm && \
cuffnorm -q  \
  --library-type fr-firststrand \
  --output-dir cuffnorm \
  --num-threads 1 \
  --labels brain,adrenal \
  cufflinks/AllSamples/merged.gtf \
  cufflinks/brain/abundances.cxb \
  cufflinks/adrenal/abundances.cxb


#-------------------------------------------------------------------------------
# STEP: differential_expression
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: differential_expression
#-------------------------------------------------------------------------------
module load mugqic/mugqic_tools/2.0.2 mugqic/R_Bioconductor/3.1.2_3.0 && \
mkdir -p DGE && \
Rscript $R_TOOLS/edger.R \
  -d design.tsv \
  -c DGE/rawCountMatrix.csv \
  -o DGE && \
Rscript $R_TOOLS/deseq.R \
  -d design.tsv \
  -c DGE/rawCountMatrix.csv \
  -o DGE


#-------------------------------------------------------------------------------
# STEP: differential_expression_goseq
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# JOB: differential_expression_goseq.dge.Contrast1
#-------------------------------------------------------------------------------
module load mugqic/mugqic_tools/2.0.2 mugqic/R_Bioconductor/3.1.2_3.0 && \
Rscript $R_TOOLS/goseq.R -p 0.1 -f 0.1 \
  -a /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Ensembl75.genes.length.tsv \
  -G /software/areas/genomics/phase2/genomes/species/Homo_sapiens.GRCh37/annotations/Homo_sapiens.GRCh37.Ensembl75.GO.tsv \
  -d DGE/Contrast1/dge_results.csv \
  -c 1,6 \
  -o DGE/Contrast1/gene_ontology_results.csv

