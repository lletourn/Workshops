# Introduction to DNA-Seq processing
This workshop will show you how to launch individual steps of a complete DNA-Seq pipeline

We will be working on a 1000 genome sample, NA12878. You can find the whole raw data on the 1000 genome website:
http://www.1000genomes.org/data

For practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

## Original Setup

The initial structure of your folders should look like this:
```
<ROOT>
|-- reference/               # genome and annotation files
|-- raw_reads/               # fastqs from the center (down sampled)
    `-- NA12878              # One sample directory
        |-- runERR_1         # Lane directory by run number. Contains the fastqs
        `-- runSRR_1         # Lane directory by run number. Contains the fastqs
`-- LabCourseNGS.csv         # sample sheet
```

### Environment setup
```
export PATH=$PATH:<TO BE DETERMINED>
```

### Software requirements

* Standalone tools:
  * [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * [BVATools](https://bitbucket.org/mugqic/bvatools/downloads)
  * [samtools](http://sourceforge.net/projects/samtools/)
  * [IGV](http://www.broadinstitute.org/software/igv/download)
  * [BWA](http://bio-bwa.sourceforge.net/)


# First data glance
So you've just received an email saying that your data is ready for download from the sequencing center of your choice. The first thing to do is download it, the second thing is making sure it is of good quality.

### Fastq files
Let's first explore the fastq file.

Try these commands
```
less -S raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz

zcat raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz | head -n4
zcat raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair2.fastq.gz | head -n4
```
From the second set of commands (the head), what was special about the output?
Why was it like that?
[Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_fastq.ex1.md)

You could also just count the reads
```
zgrep -c "^@SRR" raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz
```
Why shouldn't you just do
```
zgrep -c "^@" raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz
```
[Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/_fastq.ex2.md)


### Quality
We can't look at all the reads. Especially when working with whole genome 30x data. You could easilly have Billions of reads.

Tools like FastQC and BVATools readsqc can be used to plot many metrics from these data sets.

Let's look at the data:
```
mkdir originalQC/
java -Xmx1G -jar ~/bvatools-dev.jar readsqc --read1 raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz --read2 raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair2.fastq.gz --threads 2 --regionName SRR --output originalQC/
java -Xmx1G -jar ~/bvatools-dev.jar readsqc --read1 raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair1.fastq.gz --read2 raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair2.fastq.gz --threads 2 --regionName ERR --output originalQC/
```
Copy the images from the ```originalQC``` folder to your desktop and open the images.

What stands out in the graphs?
[Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_fastqQC.ex1.md)

All the generated graphics have their uses. This being said 2 of them are particularly useful to get an overal picture of how good or bad a run went. These are the Quality box plots and the nucleotide content graphs.

The Box plot shows the quality distribution of your data. In this case the reasons there are spikes and jumps in quality and length is because there are actually different libraries pooled together in the 2 fastq files. The sequencing lengths vary between 36,50,76 bp read lengths. The Graph goes > 100 because both ends are appended one after the other.

The quality of a base is computated using the Phread quality score.
![Phred quality score formula](../img/phredFormula.png)

The formula outputs an integer that is encoded using an [ASCII](http://en.wikipedia.org/wiki/ASCII) table. The way the lookup is done is by taking the the phred score adding 33 and using this number as a lookup in the table. The Wikipedia entry for the [FASTQ format](http://en.wikipedia.org/wiki/FASTQ_format) has a summary of the varying values.

Older illumina runs were using phred+64 instead of phred+33 to encode their fastq files.

In the SRR dataset we also see some adapters.
Why does this happen [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_fastqQC.ex2.md)


### Trimming
After this careful analysis of the raw data we see that
- Some reads have bad 3' ends.
- Some reads have adapter sequences in them.

Although nowadays this doesn't happen often, it does still happen. In some cases, miRNA, it is expected to have adapters.

Since they are not part of the genome of interest they should be removed if enough reads have them.

To be able to remove the adapters we need to feed them to a tool. In this case we will use Trimmomatic. The dapter file is already in your work folder.
We can look at the adapters
```
cat adapters.fa
```
Why are there 2 different ones? [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_trim.ex1.md)


Let's try removing them and see what happens.
```
mkdir -p reads/NA12878/runSRR_1/
mkdir -p reads/NA12878/runERR_1/

java -XX:ParallelGCThreads=1 -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair1.fastq.gz raw_reads/NA12878/runERR_1/NA12878.ERR.33.pair2.fastq.gz reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair1.fastq.gz reads/NA12878/runERR_1/NA12878.ERR.t20l32.single1.fastq.gz reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair2.fastq.gz reads/NA12878/runERR_1/NA12878.ERR.t20l32.single2.fastq.gz ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32 2> reads/NA12878/runERR_1/NA12878.ERR.trim.out

java -XX:ParallelGCThreads=1 -Xmx2G -cp $TRIMMOMATIC_JAR org.usadellab.trimmomatic.TrimmomaticPE -threads 2 -phred33 raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair1.fastq.gz raw_reads/NA12878/runSRR_1/NA12878.SRR.33.pair2.fastq.gz reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair1.fastq.gz reads/NA12878/runSRR_1/NA12878.SRR.t20l32.single1.fastq.gz reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair2.fastq.gz reads/NA12878/runSRR_1/NA12878.SRR.t20l32.single2.fastq.gz ILLUMINACLIP:adapters.fa:2:30:15 TRAILING:20 MINLEN:32 2> reads/NA12878/runSRR_1/NA12878.SRR.trim.out

cat reads/NA12878/runERR_1/NA12878.ERR.trim.out reads/NA12878/runSRR_1/NA12878.SRR.trim.out
```

What does Trimmomatic says it did? [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_trim.ex2.md)

Let's look at the graphs now

```
mkdir postTrimQC/
java -Xmx1G -jar ~/bvatools-dev.jar readsqc --read1 reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair1.fastq.gz --read2 reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair2.fastq.gz --threads 2 --regionName ERR --output postTrimQC/
java -Xmx1G -jar ~/bvatools-dev.jar readsqc --read1 reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair1.fastq.gz --read2 reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair2.fastq.gz --threads 2 --regionName SRR --output postTrimQC/
```

How does it look now? [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_trim.ex3.md)


# Alignment
The raw reads are now cleaned up of artefacts we can align each lane separatly.

Why should this be done separatly? [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_aln.ex1.md)

```
mkdir -p alignment/NA12878/runERR_1
mkdir -p alignment/NA12878/runSRR_1

bwa mem -M -t 2 -R '@RG\tID:ERR_ERR_1\tSM:NA12878\tLB:ERR\tPU:runERR_1\tCN:Broad Institute' references/b37.fasta reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair1.fastq.gz reads/NA12878/runERR_1/NA12878.ERR.t20l32.pair2.fastq.gz | java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Xmx2G -jar ${PICARD_HOME}/SortSam.jar  INPUT=/dev/stdin CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate OUTPUT=alignment/NA12878/runERR_1/NA12878.ERR.sorted.bam MAX_RECORDS_IN_RAM=500000

bwa mem -M -t 2 -R '@RG\tID:SRR_SRR_1\tSM:NA12878\tLB:SRR\tPU:runSRR_1\tCN:Broad Institute' references/b37.fasta reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair1.fastq.gz reads/NA12878/runSRR_1/NA12878.SRR.t20l32.pair2.fastq.gz | java -Djava.io.tmpdir=/lb/scratch/ -XX:ParallelGCThreads=1 -Xmx2G -jar ${PICARD_HOME}/SortSam.jar  INPUT=/dev/stdin CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate OUTPUT=alignment/NA12878/runSRR_1/NA12878.SRR.sorted.bam MAX_RECORDS_IN_RAM=500000
```

Why is it important to set Read Group information? [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_aln.ex2.md)

The details of the fields can be found in the SAM/BAM specifications [Here](http://samtools.sourceforge.net/SAM1.pdf)
For most cases, only the sample name, platform unit and library one are important. 

Why did we pipe the output of one to the other? Could we have done it differently? [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_aln.ex3.md)

We will explore the generated BAM latter.

# Lane merging
We now have alignments for each of the sequences lanes. This is not practical in it's current form. What we wan't to do now
is merge the results into one BAM.

Since we identified the reads in the BAM with read groups, even after the merging, we can still identify the origin of each read.

```
java -Xmx2G -jar ${PICARD_HOME}/MergeSamFiles.jar VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true INPUT=alignment/NA12878/runERR_1/NA12878.ERR.sorted.bam INPUT=alignment/NA12878/runSRR_1/NA12878.SRR.sorted.bam  OUTPUT=alignment/NA12878/NA12878.sorted.bam MAX_RECORDS_IN_RAM=250000 
``` 

You should now have one BAM containing all your data.
Let's double check
```
ls -l alignment/NA12878/
samtools view -H alignment/NA12878/NA12878.sorted.bam | grep "^@RG"

```

You should have your 2 read group entries.
Why did we use the ```-H``` switch? Try without. What happens? [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_merge.ex1.md)

# Cleaning up alignments
We started by cleaning up the raw reads. Now we need to fix some alignments.

The first step for this is to realign around indels and snp dense regions.
The Genome Analysis toolkit has a tool for this called IndelRealigner.

It basically runs in 2 steps
1- Find the targets
2- Realign them.

```
java -Xmx2G  -jar ${GATK_JAR} -T RealignerTargetCreator -R references/b37.fasta -o alignment/NA12878/realign.intervals -I alignment/NA12878/NA12878.sorted.bam -L 1

java -Xmx2G -jar ${GATK_JAR} -T IndelRealigner -R references/b37.fasta -targetIntervals alignment/NA12878/realign.intervals -o alignment/NA12878/NA12878.realigned.sorted.bam -I alignment/NA12878/NA12878.sorted.bam

```

How could we make this go faster? [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_realign.ex1.md)
How many regions did it think needed cleaning? [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_realign.ex2.md)

Indel Realigner also makes sure the called deletions are left aligned when there is a microsat of homopolmer.
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
one-off corrdinates and such.

This happened a lot with bwa backtrack. This happens less with bwa mem, but it still happens none the less.

```
java -Xmx2G -jar ${PICARD_HOME}/FixMateInformation.jar VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true SORT_ORDER=coordinate INPUT=alignment/NA12878/NA12878.realigned.sorted.bam OUTPUT=alignment/NA12878/NA12878.matefixed.sorted.bam MAX_RECORDS_IN_RAM=500000
```

# Mark duplicates

# Recalibration

# Extract Metrics
Once your whole bam is generated, it's always a good thing to check the data again to see if everything makes sens.

## Compute coverage
If you have data from a capture kit, you should see how well your targets worked

## Check flagstat

## Insert Size

# Variant calling

## Samtools

## GATK Unified Genotyper

## GATK Haplotyper

## Visualisation
Before jumping into IGV, we'll generate a track IGV can use to plot coverage.

Try this:

```
igvtools count -f min,max,mean alignment/NA12878/NA12878.sorted.dup.recal.bam alignment/NA12878/NA12878.sorted.dup.recal.bam.tdf b37
```

# IGV
You can get IGV [here](http://www.broadinstitute.org/software/igv/download)

Open it and choose b37 as the genome

Open your bams, the tdf we just generated should load.

What do you see...


## Aknowledgments
The format for this tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI, who I would like to thank and acknowledge. I also want to acknowledge Mathieu Bourgey, Francois Lefebvre, Maxime Caron and Guillaume Bourque for the help in building these pipelines and working with all the various datasets.
