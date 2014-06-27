# Introduction to DNA-Seq processing
This workshop will be a mix of different methods to look and explore your data.

We will be working on the same BAMs you generated from the SNV part.
Again, for practical reasons we subsampled the reads from the sample because running the whole dataset would take way too much time and resources.
This leads to some strange results in this part.

This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.

## Original Setup

The structure of your folders should now look like this:
```
<ROOT>
|-- raw_reads/                              # fastqs from the center (down sampled)
    `-- normal                              # The blood sample directory
        `-- run*_?                          # Lane directory by run number. Contains the fastqs
    `-- tumor                               # The tumor sample directory
        `-- run*_?                          # Lane directory by run number. Contains the fastqs
|-- alignment/                              # fastqs from the center (down sampled)
    `-- normal                              # The blood sample directory
        `-- normal.sorted.dup.recal.bam     # Normal alignment file
    `-- tumor                               # The tumor sample directory
        `-- tumor.sorted.dup.recal.bam      # Tumor alignment file
`-- project.nanuq.csv        # sample sheet
```

### Cheat sheets
* [Unix comand line cheat sheet](http://sites.tufts.edu/cbi/files/2013/01/linux_cheat_sheet.pdf)

### Environment setup
```
export PATH=$PATH:/home/Louis/tools/tabix-0.2.6/:/home/Louis/tools/igvtools_2.3.31/
export PICARD_HOME=/usr/local/bin
export SNPEFF_HOME=/home/Louis/tools/snpEff_v3_5_core/snpEff
export GATK_JAR=/usr/local/bin/GenomeAnalysisTK.jar
export BVATOOLS_JAR=/home/Louis/tools/bvatools-1.1/bvatools-1.1-full.jar
export TRIMMOMATIC_JAR=/usr/local/bin/trimmomatic-0.32.jar
export REF=/home/Louis/kyotoWorkshop/references/

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


# Telomeres
In this first step we will try to qualitatively see if the normal and tumor have different
telomere lengths.

One way to do this is find the telomere motif.

A good link to get various telomere repeats is the [Telomerase Database](http://telomerase.asu.edu/sequences_telomere.html)

First step, count the number of reads with these repeats.

```
```

## Aknowledgments
The format for this tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI, who I would like to thank and acknowledge. I also want to acknowledge Mathieu Bourgey, Francois Lefebvre, Maxime Caron and Guillaume Bourque for the help in building these pipelines and working with all the various datasets.
