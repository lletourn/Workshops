# Introduction
This workshop will show you how to launch individual steps of a complete RNA-Seq pipeline

The initial structure of your folders should look like this:
```
/lb/project/mugqic/hgen698_X
|-- reference/               # genome and annotation files
|-- raw_reads/               # fastqs from the center (down sampled)
    `-- SAMPLE               # One sample directory
        `-- runXXXX_Y        # Lane directory by run number. Contains the fastqs
|-- design.csv               # design file for DGE
`-- LabCourseNGS.csv         # sample sheet
```

### Be sure to set 
JOB_MAIL  == email in which to receive job completion events
WORK_DIR == Your designated work directory
```
export JOB_MAIL=dude@cool.com
export WORK_DIR=/dev/null

export MUGQIC_INSTALL_HOME=/sb/programs/analyste
export MUGQIC_INSTALL_HOME_DEV=/lb/project/mugqic/analyste_dev

export PERL5LIB=/sb/programs/analyste/software/perl5libs/lib64/perl5:${PERL5LIB}
export PERL5LIB=/sb/programs/analyste/software/perl5libs/bin/:${PERL5LIB}
export PERL5LIB=/sb/programs/analyste/software/perl5libs/share/perl5/:${PERL5LIB}
export PERL5LIB=/sb/programs/analyste/software/perl5libs/lib/perl5/:${PERL5LIB}
```

# Trimming
Why?

Look into trim.sh

Let's dig into the data with these commands:
```
mkdir -p qc
module load mugqic/bvatools/1.1

java -Xmx1G -jar $BVATOOLS_JAR readsqc --read1 raw_reads/rRNA_Dep_Brain_A/run1443_1/rRNA_Dep_Brain_A.1000003916c-B08.33.pair1.fastq.gz --read2 raw_reads/rRNA_Dep_Brain_A/run1443_1/rRNA_Dep_Brain_A.1000003916c-B08.33.pair2.fastq.gz --regionName preTrim --output qc/

java -Xmx1G -jar $BVATOOLS_JAR readsqc --read1 reads/rRNA_Dep_Brain_A/run1443_1/rRNA_Dep_Brain_A.1000003916c-B08.t30l32.pair1.fastq.gz --read2 reads/rRNA_Dep_Brain_A/run1443_1/rRNA_Dep_Brain_A.1000003916c-B08.t30l32.pair2.fastq.gz --regionName postTrim --output qc/
```

What changed?

# Alignment
Let's align

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

## Aknowledgments
The format for this tutorial has been inspired from Mar Gonzalez Porta of Embl-EBI, who I would like to thank and acknowledge.
