module load mugqic/fastqc/0.11.2
fastqc -o originalBrain raw_reads/brain/brain.pair1.fastq.gz raw_reads/brain/brain.pair2.fastq.gz
fastqc -o originalAdrenal  raw_reads/adrenal/adrenal.pair1.fastq.gz rnaSeq/raw_reads/adrenal/adrenal.pair2.fastq.gz

module load mugqic/java/openjdk-jdk1.7.0_60 mugqic/trimmomatic/0.32 && \
mkdir -p trim/brain && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 1 \
  -phred33 \
  raw_reads/brain/brain.pair1.fastq.gz \
  raw_reads/brain/brain.pair2.fastq.gz \
  trim/brain/brain.trim.pair1.fastq.gz \
  trim/brain/brain.trim.single1.fastq.gz \
  trim/brain/brain.trim.pair2.fastq.gz \
  trim/brain/brain.trim.single2.fastq.gz \
  ILLUMINACLIP:/software/areas/genomics/phase2/software/mugqic_pipeline/v1.3/lib/adapters-truseq.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/brain/brain.trim.log

#-------------------------------------------------------------------------------
# JOB: trimmomatic.adrenal
#-------------------------------------------------------------------------------
mkdir -p trim/adrenal && \
java -XX:ParallelGCThreads=1 -Xmx2G -jar $TRIMMOMATIC_JAR PE \
  -threads 1 \
  -phred33 \
  raw_reads/adrenal/adrenal.pair1.fastq.gz \
  raw_reads/adrenal/adrenal.pair2.fastq.gz \
  trim/adrenal/adrenal.trim.pair1.fastq.gz \
  trim/adrenal/adrenal.trim.single1.fastq.gz \
  trim/adrenal/adrenal.trim.pair2.fastq.gz \
  trim/adrenal/adrenal.trim.single2.fastq.gz \
  ILLUMINACLIP:/software/areas/genomics/phase2/software/mugqic_pipeline/v1.3/lib/adapters-truseq.fa:2:30:15 \
  TRAILING:30 \
  MINLEN:32 \
  2> trim/adrenal/adrenal.trim.log

module load mugqic/fastqc/0.11.2
fastqc -o trimmedBrain trim/brain/brain.trim.pair1.fastq.gz trim/brain/brain.trim.pair2.fastq.gz
fastqc -o trimmedAdrenal trim/adrenal/adrenal.trim.pair1.fastq.gz trim/adrenal/adrenal.trim.pair2.fastq.gz

#-------------------------------------------------------------------------------
# STEP: tophat
#-------------------------------------------------------------------------------
echo "module load mugqic/bowtie/2.1.0 mugqic/tophat/2.0.10 mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
mkdir -p alignment/brain/tophat/ && \
tophat --no-coverage-search \
 --rg-library \"brain\" \
 --rg-platform \"ILLUMINA\" \
 --rg-platform-unit \"2\" \
 --rg-sample \"brain\" \
 --rg-id brain \
 --library-type fr-firststrand \
 -o alignment/brain/tophat/ -p 2 \
 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/genome/bowtie2_index/Homo_sapiens.GRCh37 \
 trim/brain/brain.trim.pair1.fastq.gz \
 trim/brain/brain.trim.pair2.fastq.gz && \
 java -XX:ParallelGCThreads=1 -Xmx2G -jar \${PICARD_HOME}/MergeSamFiles.jar \
 VALIDATION_STRINGENCY=SILENT \
 CREATE_INDEX=true \
 INPUT=alignment/brain/tophat/accepted_hits.bam \
 INPUT=alignment/brain/tophat/unmapped.bam \
 OUTPUT=alignment/brain/tophat/brain.sorted.bam \
 SORT_ORDER=coordinate" | qsub -V -m ae -M $JOB_MAIL -d . -j oe -N brain -l walltime=24:00:00 -q metaq -l nodes=1:ppn=2

#-------------------------------------------------------------------------------
# STEP: tophat
#-------------------------------------------------------------------------------

echo "module load mugqic/bowtie/2.1.0 mugqic/tophat/2.0.10 mugqic/samtools/0.1.19 mugqic/java/openjdk-jdk1.7.0_60 mugqic/picard/1.123 && \
mkdir -p alignment/adrenal/tophat/ && \
tophat --no-coverage-search \
 --rg-library \"adrenal\" \
 --rg-platform \"ILLUMINA\" \
 --rg-platform-unit \"2\" \
 --rg-sample \"adrenal\" \
 --rg-id adrenal \
 --library-type fr-firststrand \
 -o alignment/adrenal/tophat/ -p 2 \
 $MUGQIC_INSTALL_HOME/genomes/species/Homo_sapiens.GRCh37/genome/bowtie2_index/Homo_sapiens.GRCh37 \
 trim/adrenal/adrenal.trim.pair1.fastq.gz \
 trim/adrenal/adrenal.trim.pair2.fastq.gz && \
 java -XX:ParallelGCThreads=1 -Xmx2G -jar ${PICARD_HOME}/MergeSamFiles.jar \
 VALIDATION_STRINGENCY=SILENT \
 CREATE_INDEX=true \
 INPUT=alignment/adrenal/tophat/accepted_hits.bam \
 INPUT=alignment/adrenal/tophat/unmapped.bam \
 OUTPUT=alignment/adrenal/tophat/adrenal.sorted.bam \
 SORT_ORDER=coordinate"  | qsub -V -m ae -M $JOB_MAIL -d . -j oe -N adrenal -l walltime=24:00:00 -q metaq -l nodes=1:ppn=2

