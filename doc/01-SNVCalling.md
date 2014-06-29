# Introduction to DNA-Seq processing
This workshop will show you how to launch individual steps of a complete DNA-Seq pipeline for Cancer analysis

We will be working on a CageKid sample pair, patient C0098.
The CageKid project is part of ICGC and is focused on renal cancer in many of it's forms.
The raw data can be found on EGA and calls, RNA and DNA, can be found on the ICGC portal.
More details about CageKid here:
http://www.cng.fr/cagekid/

For practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

## Original Setup

The initial structure of your folders should look like this:
```
<ROOT>
|-- raw_reads/               # fastqs from the center (down sampled)
    `-- normal               # The blood sample directory
        `-- run*_?           # Lane directory by run number. Contains the fastqs
    `-- tumor                # The tumor sample directory
        `-- run*_?           # Lane directory by run number. Contains the fastqs
`-- project.nanuq.csv        # sample sheet
```

### Cheat sheets
* [Unix comand line cheat sheet](http://sites.tufts.edu/cbi/files/2013/01/linux_cheat_sheet.pdf)

### Environment setup
```
export APP_ROOT=/home/training/Applications/
export PATH=$PATH:$APP_ROOT/bwa-0.7.9a:$APP_ROOT/tabix-0.2.6/:$APP_ROOT/IGVTools
export PICARD_HOME=$APP_ROOT/picard-tools-1.115/
export SNPEFF_HOME=$APP_ROOT/snpEff/
export GATK_JAR=$APP_ROOT/gatk/GenomeAnalysisTK.jar
export BVATOOLS_JAR=$APP_ROOT/bvatools-1.3/bvatools-1.3-full.jar
export TRIMMOMATIC_JAR=$APP_ROOT/Trimmomatic-0.32/trimmomatic-0.32.jar
export STRELKA_HOME=$APP_ROOT/strelka-1.0.13/
export REF=/home/training/ebiCancerWorkshop201407/references/

cd $HOME/ebiCancerWorkshop201407
```

### Software requirements
These are all already installed, but here are the original links.

  * [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [BVATools](https://bitbucket.org/mugqic/bvatools/downloads)
  * [SAMTools](http://sourceforge.net/projects/samtools/)
  * [IGV](http://www.broadinstitute.org/software/igv/download)
  * [BWA](http://bio-bwa.sourceforge.net/)
  * [Genome Analysis Toolkit](http://www.broadinstitute.org/gatk/)
  * [Picard](http://picard.sourceforge.net/)
  * [SnpEff](http://snpeff.sourceforge.net/)
  * [MuTect](http://www.broadinstitute.org/cancer/cga/mutect)
  * [Strelka](https://sites.google.com/site/strelkasomaticvariantcaller/)


# First data glance
So you've just received an email saying that your data is ready for download from the sequencing center of your choice.
The first thing to do is download it, the second thing is making sure it is of good quality.

### Fastq files
Let's first explore the fastq file.

Try these commands
```
zless -S raw_reads/normal/runD0YR4ACXX_1/normal.64.pair1.fastq.gz

zcat raw_reads/normal/runD0YR4ACXX_1/normal.64.pair1.fastq.gz | head -n4
zcat raw_reads/normal/runD0YR4ACXX_1/normal.64.pair2.fastq.gz | head -n4
```
From the second set of commands (the head), what was special about the output?
Why was it like that?
[Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_fastq.ex1.md)

You could also just count the reads
```
zgrep -c "^@HISEQ2" raw_reads/normal/runD0YR4ACXX_1/normal.64.pair1.fastq.gz
```
Why shouldn't you just do
```
zgrep -c "^@" raw_reads/normal/runD0YR4ACXX_1/normal.64.pair1.fastq.gz
```
[Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_fastq.ex2.md)


### Quality
We can't look at all the reads. Especially when working with whole genome 50x data. You could easilly have Billions of reads.

Tools like FastQC and BVATools readsqc can be used to plot many metrics from these data sets.

Let's look at the data:

```
# Generate original QC
mkdir originalQC/
java7 -Xmx1G -jar ${BVATOOLS_JAR} readsqc --quality 64 \
  --read1 raw_reads/normal/runD0YR4ACXX_1/normal.64.pair1.fastq.gz \
  --read2 raw_reads/normal/runD0YR4ACXX_1/normal.64.pair2.fastq.gz \
  --threads 2 --regionName normalD0YR4ACXX_1 --output originalQC/
```

Open the images.

What stands out in the graphs?
[Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_fastqQC.ex1.md)

All the generated graphics have their uses. This being said 2 of them are particularly useful to get an overal picture of how good or bad a run went. These are the Quality box plots and the nucleotide content graphs.

The quality of a base is computated using the Phread quality score.
![Phred quality score formula](../img/phredFormula.png)

The formula outputs an integer that is encoded using an [ASCII](http://en.wikipedia.org/wiki/ASCII) table. The way the lookup is done is by taking the the phred score adding 33 and using this number as a lookup in the table. The Wikipedia entry for the [FASTQ format](http://en.wikipedia.org/wiki/FASTQ_format) has a summary of the varying values.

Older illumina runs, and the data here, used phred+64 instead of phred+33 to encode their fastq files.

We see a little bit of adapter in the sequences.
Why does this happen [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_fastqQC.ex2.md)


### Trimming
After this careful analysis of the raw data we see that
- Some reads have bad 3' ends.
- Some reads have adapter sequences in them.
- Data needs to be converted into phred+33 from phred+64

Although nowadays this doesn't happen often, it does still happen. In some cases, miRNA, it is expected to have adapters.

Since they are not part of the genome of interest they should be removed if enough reads have them.

To be able to remove the adapters we need to feed them to a tool. In this case we will use Trimmomatic. The adapter file is already in your work folder.
We can look at the adapters
```
cat adapters.fa
```
Why are there 2 different ones? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_trim.ex1.md)

Another reason we want to run Trimmomatic here is to convert the data from phred+33 to phred+64. Trimmomatic does this directly.
In modern datasets this is not needed.


Let's try removing them and see what happens.
```
# Trim and convert data
for file in raw_reads/*/run*_?/*.pair1.fastq.gz;
do
  FNAME=`basename $file`;
  DIR=`dirname $file`;
  OUTPUT_DIR=`echo $DIR | sed 's/raw_reads/reads/g'`;

  mkdir -p $OUTPUT_DIR;
  java7 -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred64 \
    $file \
    ${file%.pair1.fastq.gz}.pair2.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.pair1.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.single1.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.pair2.fastq.gz \
    ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.t30l50.single2.fastq.gz \
    TOPHRED33 ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:30 MINLEN:50 \
    2> ${OUTPUT_DIR}/${FNAME%.64.pair1.fastq.gz}.trim.out ; 
done

cat reads/normal/runD0YR4ACXX_1/normal.trim.out
```

What does Trimmomatic says it did? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_trim.ex2.md)

Since the data was so good to start with, we won't regenerate the graph post trim (or you could do it as an extra exercise if you wish)

# Alignment
The raw reads are now cleaned up of artefacts we can align each lane separatly.

Why should this be done separatly? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_aln.ex1.md)

```
# Align data
for file in reads/*/run*_?/*.pair1.fastq.gz;
do
  FNAME=`basename $file`;
  DIR=`dirname $file`;
  OUTPUT_DIR=`echo $DIR | sed 's/reads/alignment/g'`;
  SNAME=`echo $file | sed 's/reads\/\([^/]\+\)\/.*/\1/g'`;
  RUNID=`echo $file | sed 's/.*\/run\([^_]\+\)_.*/\1/g'`;
  LANE=`echo $file | sed 's/.*\/run[^_]\+_\(.\).*/\1/g'`;

  mkdir -p $OUTPUT_DIR;

  bwa mem -M -t 3 \
    -R "@RG\\tID:${SNAME}_${RUNID}_${LANE}\\tSM:${SNAME}\\tLB:${SNAME}\\tPU:${RUNID}_${LANE}\\tCN:Centre National de Genotypage\\tPL:ILLUMINA" \
    ${REF}/bwa/b37.fasta \
    $file \
    ${file%.pair1.fastq.gz}.pair2.fastq.gz \
  | java7 -Xmx2G -jar ${PICARD_HOME}/SortSam.jar \
    INPUT=/dev/stdin \
    OUTPUT=${OUTPUT_DIR}/${SNAME}.sorted.bam \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000
done
```

Why is it important to set Read Group information? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_aln.ex2.md)

The details of the fields can be found in the SAM/BAM specifications [Here](http://samtools.sourceforge.net/SAM1.pdf)
For most cases, only the sample name, platform unit and library one are important. 

Why did we pipe the output of one to the other? Could we have done it differently? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_aln.ex3.md)

We will explore the generated BAM latter.

# Lane merging
We now have alignments for each of the sequences lanes. This is not practical in it's current form. What we wan't to do now
is merge the results into one BAM.

Since we identified the reads in the BAM with read groups, even after the merging, we can still identify the origin of each read.

```
# Merge Data
java7 -Xmx2G -jar ${PICARD_HOME}/MergeSamFiles.jar \
  INPUT=alignment/normal/runC0LWRACXX_1/normal.sorted.bam \
  INPUT=alignment/normal/runC0LWRACXX_6/normal.sorted.bam \
  INPUT=alignment/normal/runC0PTAACXX_6/normal.sorted.bam \
  INPUT=alignment/normal/runC0PTAACXX_7/normal.sorted.bam \
  INPUT=alignment/normal/runC0PTAACXX_8/normal.sorted.bam \
  INPUT=alignment/normal/runC0R2BACXX_6/normal.sorted.bam \
  INPUT=alignment/normal/runC0R2BACXX_7/normal.sorted.bam \
  INPUT=alignment/normal/runC0R2BACXX_8/normal.sorted.bam \
  INPUT=alignment/normal/runD0YR4ACXX_1/normal.sorted.bam \
  INPUT=alignment/normal/runD0YR4ACXX_2/normal.sorted.bam \
  OUTPUT=alignment/normal/normal.sorted.bam \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

java7 -Xmx2G -jar ${PICARD_HOME}/MergeSamFiles.jar \
  INPUT=alignment/tumor/runBC0TV0ACXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0LVJACXX_6/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0PK4ACXX_7/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0PK4ACXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0R29ACXX_7/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0R29ACXX_8/tumor.sorted.bam \
  INPUT=alignment/tumor/runC0TTBACXX_3/tumor.sorted.bam \
  INPUT=alignment/tumor/runD114WACXX_8/tumor.sorted.bam \
  OUTPUT=alignment/tumor/tumor.sorted.bam \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

``` 

You should now have one normal and one tumor BAM containing all your data.
Let's double check
```
ls -l alignment/normal/
samtools view -H alignment/normal/normal.sorted.bam | grep "^@RG"

```

You should have your 10 read group entries.
Why did we use the ```-H``` switch? Try without. What happens? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_merge.ex1.md)

## SAM/BAM
Let's spend some time to explore bam files.

try
```
samtools view alignment/normal/normal.sorted.bam | head -n2
```

Here you have examples of alignment results.
A full description of the flags can be found in the SAM specification
http://samtools.sourceforge.net/SAM1.pdf

You can try using picards explain flag site to understand what is going on with your reads
http://picard.sourceforge.net/explain-flags.html

The flag is the 2nd column.


You can use samtools to filter.

```
# Say you want to count the *un-aligned* reads you can use
samtools view -c -f4 alignment/normal/normal.sorted.bam

# Or you want to count the *aligned* reads you can use
samtools view -c -F4 alignment/normal/normal.sorted.bam
```

We won't go into too much detail at this point since we want to concentrate on cancer specific issues now.

How many reads mapped and unmapped were there? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_sambam.ex3.md)


Another useful bit of information in the SAM is the CIGAR string.
It's the 6th column in the file. This column explains how the alignment was achieved.
M == base aligns *but doesn't have to be a match*. A SNP will have an M even if it disagrees with the reference.
I == Insertion
D == Deletion
S == soft-clips. These are handy to find un removed adapters, viral insertions, etc.

An in depth explanation of the CIGAR can be found [here](http://genome.sph.umich.edu/wiki/SAM)
The exact details of the cigar string can be found in the SAM spec as well.
Another good site

# Cleaning up alignments
We started by cleaning up the raw reads. Now we need to fix some alignments.

The first step for this is to realign around indels and snp dense regions.
The Genome Analysis toolkit has a tool for this called IndelRealigner.

It basically runs in 2 steps
1- Find the targets
2- Realign them.


For cancer there is a subtility

```
# Realign
java7 -Xmx2G  -jar ${GATK_JAR} \
  -T RealignerTargetCreator \
  -R ${REF}/b37.fasta \
  -o alignment/normal/realign.intervals \
  -I alignment/normal/normal.sorted.bam \
  -I alignment/tumor/tumor.sorted.bam \
  -L 19

java7 -Xmx2G -jar ${GATK_JAR} \
  -T IndelRealigner \
  -R ${REF}/b37.fasta \
  -targetIntervals alignment/normal/realign.intervals \
  --nWayOut .realigned.bam \
  -I alignment/normal/normal.sorted.bam \
  -I alignment/tumor/tumor.sorted.bam

  mv normal.sorted.realigned.bam alignment/normal/
  mv tumor.sorted.realigned.bam alignment/tumor/

```

Why did we use both normal and tumor together? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_realign.ex3.md)

How could we make this go faster? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_realign.ex1.md)
How many regions did it think needed cleaning? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_realign.ex2.md)

Indel Realigner also makes sure the called deletions are left aligned when there is a microsatellite or homopolymer.
```
This
ATCGAAAA-TCG
into
ATCG-AAAATCG

or
ATCGATATATATA--TCG
into
ATCG--ATATATATATCG
```

This makes it easier for down stream tools.

# FixMates
This step shouldn't be necessary...but it is.

This goes through the BAM file and find entries which don't have their mate information written properly.

This used to be a problem in the GATKs realigner, but they fixed it. It shouldn't be a problem with aligners like BWA, but there are always corner cases that create
one-off coordinates and such.

This happened a lot with bwa backtrack. This happens less with bwa mem, but it still happens none the less.

```
# Fix Mate
java7 -Xmx2G -jar ${PICARD_HOME}/FixMateInformation.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000 \
  INPUT=alignment/normal/normal.sorted.realigned.bam \
  OUTPUT=alignment/normal/normal.matefixed.bam
java7 -Xmx2G -jar ${PICARD_HOME}/FixMateInformation.jar \
  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=500000 \
  INPUT=alignment/tumor/tumor.sorted.realigned.bam \
  OUTPUT=alignment/tumor/tumor.matefixed.bam
```

# Mark duplicates
As the step says, this is to mark duplicate reads.
What are duplicate reads? What are they caused by? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_markdup.ex1.md)

What are the ways to detect them? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_markdup.ex2.md)

Here we will use picards approach:

```
# Mark Dups
java7 -Xmx2G -jar ${PICARD_HOME}/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/normal/normal.matefixed.bam \
  OUTPUT=alignment/normal/normal.sorted.dup.bam \
  METRICS_FILE=alignment/normal/normal.sorted.dup.metrics

java7 -Xmx2G -jar ${PICARD_HOME}/MarkDuplicates.jar \
  REMOVE_DUPLICATES=false CREATE_MD5_FILE=true VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
  INPUT=alignment/tumor/tumor.matefixed.bam \
  OUTPUT=alignment/tumor/tumor.sorted.dup.bam \
  METRICS_FILE=alignment/tumor/tumor.sorted.dup.metrics
```

We can look in the metrics output to see what happened.
We can see that it computed seperate measures for each library.
Why is this important to do and not combine everything? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_markdup.ex3.md)

How many duplicates were there? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_markdup.ex4.md)

This is pretty spot on for a whole genome project. <5% should be expected.

# Recalibration
This is the last BAM cleaning up step.

The goal for this step is to try to recalibrate base quality scores. The vendors tend to inflate the values of the bases in the reads.
Also, this step tries to lower the scores of some biased motifs for some technologies.

It runs in 2 steps, 
1- Build covariates based on context and known snp sites
2- Correct the reads based on these metrics

```
# Recalibrate
for i in normal tumor
do
  java7 -Xmx2G -jar ${GATK_JAR} \
    -T BaseRecalibrator \
    -nct 2 \
    -R ${REF}/b37.fasta \
    -knownSites ${REF}/dbSnp-137.vcf.gz \
    -L 19:50500000-52502000 \
    -o alignment/${i}/${i}.sorted.dup.recalibration_report.grp \
    -I alignment/${i}/${i}.sorted.dup.bam

    java7 -Xmx2G -jar ${GATK_JAR} \
      -T PrintReads \
      -nct 2 \
      -R ${REF}/b37.fasta \
      -BQSR alignment/${i}/${i}.sorted.dup.recalibration_report.grp \
      -o alignment/${i}/${i}.sorted.dup.recal.bam \
      -I alignment/${i}/${i}.sorted.dup.bam
done
```

Just to see how things change let's make GATK recalibrate after a first pass
```
# Check Recalibration
for i in normal tumor
do
  java7 -Xmx2G -jar ${GATK_JAR} \
    -T BaseRecalibrator \
    -nct 2 \
    -R ${REF}/b37.fasta \
    -knownSites ${REF}/dbSnp-137.vcf.gz \
    -L 19:50500000-52502000 \
    -o alignment/${i}/${i}.sorted.dup.recalibration_report.seconnd.grp \
    -I alignment/${i}/${i}.sorted.dup.bam \
    -BQSR alignment/${i}/${i}.sorted.dup.recalibration_report.grp

  java7 -Xmx2G -jar ${GATK_JAR} \
    -T AnalyzeCovariates \
    -R ${REF}/b37.fasta \
    -before alignment/${i}/${i}.sorted.dup.recalibration_report.grp \
    -after alignment/${i}/${i}.sorted.dup.recalibration_report.seconnd.grp \
    -csv alignment/${i}/BQSR.${i}.csv \
    -plots alignment/${i}/BQSR.${i}.pdf
done
```

The graphs don't mean much because we downsampled the data quite a bit. With a true whole genome or whole exome dataset we can see a bigger effect.

# Extract Metrics
Once your whole bam is generated, it's always a good thing to check the data again to see if everything makes sens.

## Compute coverage
If you have data from a capture kit, you should see how well your targets worked

Both GATK and BVATools have depth of coverage tools. We wrote our own in BVAtools because
- GATK was deprecating theirs, but they changed their mind
- GATK's is very slow
- We were missing some output that we wanted from the GATK's one (GC per interval, valid pairs, etc)

Here we'll use the GATK one

```
# Get Depth
for i in normal tumor
do
  java7  -Xmx2G -jar ${GATK_JAR} \
    -T DepthOfCoverage \
    --omitDepthOutputAtEachBase \
    --summaryCoverageThreshold 10 \
    --summaryCoverageThreshold 25 \
    --summaryCoverageThreshold 50 \
    --summaryCoverageThreshold 100 \
    --start 1 --stop 500 --nBins 499 -dt NONE \
    -R ${REF}/b37.fasta \
    -o alignment/${i}/${i}.sorted.dup.recal.coverage \
    -I alignment/${i}/${i}.sorted.dup.recal.bam \
    -L 19:50500000-52502000 &
done
wait

# Look at the coverage
less -S alignment/normal/normal.sorted.dup.recal.coverage.sample_interval_summary
less -S alignment/tumor/tumor.sorted.dup.recal.coverage.sample_interval_summary
```

Coverage is the expected ~70-110x. 
summaryCoverageThreshold is a usefull function to see if your coverage is uniform.
Another way is to compare the mean to the median. If both are almost equal, your coverage is pretty flat. If both are quite different
That means something is wrong in your coverage. A mix of WGS and WES would show very different mean and median values.

## Insert Size
```
# Get insert size
for i in normal tumor
do
  java -Xmx2G -jar ${PICARD_HOME}/CollectInsertSizeMetrics.jar \
    VALIDATION_STRINGENCY=SILENT \
    REFERENCE_SEQUENCE=${REF}/b37.fasta \
    INPUT=alignment/${i}/${i}.sorted.dup.recal.bam \
    OUTPUT=alignment/${i}/${i}.sorted.dup.recal.metric.insertSize.tsv \
    HISTOGRAM_FILE=alignment/${i}/${i}.sorted.dup.recal.metric.insertSize.histo.pdf \
    METRIC_ACCUMULATION_LEVEL=LIBRARY
done

#look at the output
less -S alignment/normal/normal.sorted.dup.recal.metric.insertSize.tsv
less -S alignment/tumor/tumor.sorted.dup.recal.metric.insertSize.tsv
```

## Alignment metrics
For the alignment metrics, we used to use ```samtools flagstat``` but with bwa mem since some reads get broken into pieces, the numbers are a bit confusing.
You can try it if you want.

We prefer the Picard way of computing metrics

```
# Get alignment metrics
for i in normal tumor
do
  java -Xmx2G -jar ${PICARD_HOME}/CollectAlignmentSummaryMetrics.jar \
    VALIDATION_STRINGENCY=SILENT \
    REFERENCE_SEQUENCE=${REF}/b37.fasta \
    INPUT=alignment/${i}/${i}.sorted.dup.recal.bam \
    OUTPUT=alignment/${i}/${i}.sorted.dup.recal.metric.alignment.tsv \
    METRIC_ACCUMULATION_LEVEL=LIBRARY
done

# explore the results
less -S alignment/normal/normal.sorted.dup.recal.metric.alignment.tsv
less -S alignment/tumor/tumor.sorted.dup.recal.metric.alignment.tsv

```


# Variant calling
Here we will try 3 variant callers.
- SAMtools
- MuTecT
- Strelka

Other candidates
- Varscan 2
- Virmid
- Somatic sniper

many, MANY others can be found here:
https://www.biostars.org/p/19104/

In our case, let's start with:
```
mkdir pairedVariants
```

## SAMtools
```
# Variants SAMTools
samtools mpileup -L 1000 -B -q 1 -D -S -g \
  -f ${REF}/b37.fasta \
  -r 19:50500000-52502000 \
  alignment/normal/normal.sorted.dup.recal.bam \
  alignment/tumor/tumor.sorted.dup.recal.bam \
  | bcftools view -vcg -T pair - \
  > pairedVariants/mpileup.vcf
```

## Broad MuTecT
```
# Variants MuTecT
# Note MuTecT only works with Java 6, 7 will give you an error
# if you get "Comparison method violates its general contract!
# you used java 7"
java -Xmx2G -jar ${MUTECT_JAR} \
  -T MuTect \
  -R ${REF}/b37.fasta \
  -dt NONE -baq OFF --validation_strictness LENIENT -nt 2 \
  --dbsnp ${REF}/dbSnp-137.vcf \
  --cosmic ${REF}/b37_cosmic_v54_120711.vcf \
  --input_file:normal alignment/normal/normal.sorted.dup.recal.bam \
  --input_file:tumor alignment/tumor/tumor.sorted.dup.recal.bam \
  --out pairedVariants/mutect.call_stats.txt \
  --coverage_file pairedVariants/mutect.wig.txt \
  -pow pairedVariants/mutect.power \
  -vcf pairedVariants/mutect.vcf \
  -L 19:50500000-52502000
```

## Illumina Strelka
```
# Variants Strelka
cp ${STRELKA_HOME}/etc/strelka_config_bwa_default.ini ./
# Fix ini since we subsampled
sed 's/isSkipDepthFilters =.*/isSkipDepthFilters = 1/g' -i strelka_config_bwa_default.ini

${STRELKA_HOME}/bin/configureStrelkaWorkflow.pl \
  --normal=alignment/normal/normal.sorted.dup.recal.bam \
  --tumor=alignment/tumor/tumor.sorted.dup.recal.bam \
  --ref=${REF}/b37.fasta \
  --config=strelka_config_bwa_default.ini \
  --output-dir=pairedVariants/strelka/

  cd pairedVariants/strelka/
  make -j3
  cd $HOME/ebiCancerWorkshop201407

  cp pairedVariants/strelka/results/passed.somatic.snvs.vcf pairedVariants/strelka.vcf
```

Now we have variants from all three methods. Let's compress and index the vcfs for futur visualisation.
```
for i in pairedVariants/*.vcf;do bgzip -c $i > $i.gz ; tabix -p vcf $i.gz;done
```

Let's look at a compressed vcf.
```
zless -S variants/mpileup.vcf.gz
```

Details on the spec can be found here:
http://vcftools.sourceforge.net/specs.html

Fields vary from caller to caller. Some values are more constant.
The ref vs alt alleles, variant quality (QUAL column) and the per-sample genotype (GT) values are almost always there.

# Annotations
We typically use snpEff but many use annovar and VEP as well.

Let's run snpEff
```
# SnpEff
java7  -Xmx6G -jar ${SNPEFF_HOME}/snpEff.jar \
  eff -v -c ${SNPEFF_HOME}/snpEff.config \
  -o vcf \
  -i vcf \
  -stats pairedVariants/mpileup.snpeff.vcf.stats.html \
  hg19 \
  pairedVariants/mpileup.vcf \
  > pairedVariants/mpileup.snpeff.vcf

less -S pairedVariants/mpileup.snpeff.vcf
```
We can see in the vcf that snpEff added a few sections. These are hard to decipher directly from the VCF other tools or scripts,
need to be used to make sens of this.

For now we will skip this step since you will be working with gene annotations in your next workshop.

Take a look at the HTML stats file snpEff created. It contains some metrics on the variants it analysed.

## Visualisation
Before jumping into IGV, we'll generate a track IGV can use to plot coverage.

Try this:

```
# Coverage Track
for i in normal tumor
do
  igvtools count \
    -f min,max,mean \
    alignment/${i}/${i}.sorted.dup.recal.bam \
    alignment/${i}/${i}.sorted.dup.recal.bam.tdf \
    b37
done
```

# IGV
You can get IGV [here](http://www.broadinstitute.org/software/igv/download)

Open it and choose b37 as the genome

Open your BAM file, the tdf we just generated should load.
Load your vcfs as well.

Find an indel. What's different between the snp callers? [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_vis.ex1.md)
Go to 19:50500000-52502000 what is interesting here?  [Solution](https://github.com/lletourn/Workshops/blob/ebiCancerWorkshop201407/solutions/_vis.ex2.md)

Look around...


## Aknowledgments
The format for this tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI, who I would like to thank and acknowledge. I also want to acknowledge Mathieu Bourgey, Francois Lefebvre, Maxime Caron and Guillaume Bourque for the help in building these pipelines and working with all the various datasets.
