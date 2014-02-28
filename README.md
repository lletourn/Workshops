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
The data is now cleaned up of artefacts we can align each lane seperatly.

Why should this be done seperatly? [Solution](https://github.com/lletourn/Workshops/blob/kyoto201403/blob/solutions/_aln.ex1.md)

```
```

## TopHat
Look into align.sh

If we would have aligned to a transcriptome what would be different?

## MarkDup
Look into mdup.sh
Run all the lines but the last.

Why markdup?

Should we?

Look in the metrics file. Are there any dups?
Does it fit with the QC?

Let's try the last command we skipped from mdup.sh

## Visualisation
We'll need a coverage marker to find our way.

Try this:

```
module load mugqic/igvtools/2.3.14
igvtools count -f min,max,mean alignment/rRNA_Dep_Brain_A/rRNA_Dep_Brain_A.merged.mdup.bam alignment/rRNA_Dep_Brain_A/rRNA_Dep_Brain_A.merged.mdup.bam.tdf b37
```

### IGV
You can get IGV [here](http://www.broadinstitute.org/software/igv/download)

Open it and choose b37 as the genome

Open your bams, the tdf should load.

What do you see...

## RNA QC
Ok the command should be done now.

## Counting reads
IGV is not meant for this, but this is needed to understand DGE.

We'll use [*htseq*](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)

But it requires special ordering.

Look at count.sh

run the first lines until htseq.

Why did we need to run these sorts?

Now run the last lines, the htseq commands.

## Saturation
Why is this useful?

Look and run saturation.sh

What do we see in the images?

## Cufflinks

## DGE
A good paper that explains the various methods can be found here:
Dillies, M.-A., Rau, A., Aubert, J., Hennequet-Antier, C., Jeanmougin, M., Servant, N., Keime, C., et al. (2013). A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis. Briefings in bioinformatics, 14(6), 671–83. doi:10.1093/bib/bbs046

Run edger
Run DESeq



## Aknowledgments
The format for this tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI, who I would like to thank and acknowledge. I also want to acknowledge Mathieu Bourgey, Francois Lefebvre, Maxime Caron and Guillaume Bourque for the help in building these pipelines and working with all the various datasets.
