# Where's my data!
You should find your data on the USB key I'll pass around.

One 1000genome BAM, one region.bed

The first thing to do is load it.

- Login to galaxy
- On the left pane you have a 'Get Data' header, click it
- Choose 'Upload File'

We need to setup information for Galaxy to process the file correctly.
We will be uploaded a 1000genome sub-sampled dataset in BAM format.

Choose BAM format and upload your BAM 'NA12878.sorted.bam'

pick the right Genome, ours is Human b37.
What's the difference between these?


# BAM? I want fastqs
The problem is, many sites, like ours, gives the data in BAM (Binary Alignment Map) format.
Since we want to try out the steps manually, we need to transform the data from BAM to Fastqs.

- Under NGS Tools choose Picard->Sam To Fastq
- Make sure your BAM is selected and that the checkbox for paired reads is selected as well

# QC, always QC
You need to QC your reads.
Your facility might offer QC graphs, we do, but they might not. Here let's look at it using FastQC.

- Under NGS Tools, pick QC and manipulation, choose FastQC Comprehensive
- Run it on *both* fastqs to compare

What do we see?

# Alignment
Data seems ok, let's align it.

There are abunch of aligners out there. Most were developped for a specific functions, contig alignment, long read alignment, high read error alignment, etc.
For NGS, most groups fall back on BWA to accomplish this task.

Here we will do the same.

Reads in this dataset are short, so we should use BWA backtrack/aln, but typically mem is used with reads >= 100bp. We will use mem.

- Under NGS Tools, pick NGS Mapping
- Pick BWA for Illumina Reads

Why are there 2 BWA choices?

Try to fill out the different fields.
Here are some hints:
- our data is paired end
- We want to align to hg19

Don't pick 'Common Parameters', pick full list.

We will keep the parameters as-is, except for Read Group.

*ALWAYS* set the read group. It it optional but for debugging purposes later it is always good practice to track where the data came from.
Also, some tools, like the Broads Genome Analysis Toolkit, will refuse to work with your BAMs if they don't have the Read Group Set.

More details on the Read Group and the SAM format in general here:
http://samtools.github.io/hts-specs/

Checkout SAMv1

for the purposes of the course you can put anything in the read group fields, but usually you would set the Sample, Library and Platform Unit fields at a minimum.

## Alignment data
Click on the 'eye' (view data)

this is what a SAM (readable BAM) looks like. Again details on the format can be found in the spec.

Let's look at a few flag settings.
This site https://broadinstitute.github.io/picard/explain-flags.html is very useful to help understand what these mean.

## SAM-to-BAM and sort
SAM is nice to look at but it's not practical for it's size and access times.

Typically we give the ouput of BWA directly to Picard or Samtools to generate a BAM, sort it and index it in one shot.

Here we'll do it step by step.

- Go into NGS Tools, Pick NGS Sam Tools SAM-to-BAM
- Pick your SAM and launch the tool

Here the wrapper does the extra work of sorting the BAM at the same time. So you are left with a Coordinate sorted BAM.

What other sorting methods could be useful and why?

## Mark Dups
We now need to mark duplicates.

Why should we do this?

- Under NGS Tools, pick NGS Picard
- Mark Duplicates
- Un check Remove duplicates from output file

What will happen to dups since we don't remove them? Why don't we remove them?

You can look at the output metrics to see how many dups there were.

# Alignment stats
First off, in a normal processing pipeline I would run many other steps that I will explain after all this. But for the sake of time, let's pretend we ran everything and
we have a pristine BAM

To generate most stats we will be using the Picard suite of tools. Most outputs are in text form but there are a few graphs for some of the stats (insert size for example)

- Under NGS Tools, pick NGS Picard
- Run Alignment summary
- Run Insert Size distribution
- Under NGS SAM Tools choose Analyse SAM/BAM, choose region, choose your region.bed file.

# Coverage
We also want to see how uniform our coverage was. This is particularly useful with Exon capture kits (or any other capture kit).

- Under NGS Tools, pick NGS GATK (beta)
- Depth of coverage
- under Summary coverage thresholds add 5x and 10x by clicking the 'Add button'
- Also pick advanced parameters and under -L (interval list) add the region.bed file.
- We will also generate a bigWig under DeepTools choose bamCoverage

# Calling variants
Of course with DNASeq what you want is to call variants. We will try with 2 callers here, Freebayes and samtools

- Under NGS Tools, pick NGS SAM tools, MPileup, generate likely hoods
- Once this is done go back to SAM tolls and choose bcftools view

Now we have a VCF to explore

# Visualization
We will now explore this data in IGV in person.

