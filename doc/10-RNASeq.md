# Introduction
This workshop will show you how to launch individual steps of a complete RNA-Seq pipeline

# Copy the data
We moved all the RNA data in the 'class10' folder.

```
cp -r /home/class10/rnaseq ./
```

The initial structure of your folders should now look like this:
```
/home/classXX
|-- raw_reads/                        # fastqs from the center (down sampled)
    `-- brain                         # One sample directory
        `-- brain.pair[12].fastq.gz   # Fastqs of the sample
    `-- adrenal                       # One sample directory
        `-- adrenal.pair[12].fastq.gz # Fastqs of the sample
`-- design.csv                        # design file for DGE
`-- steps.batch.1.sh                  # First part
`-- steps.batch.2.sh                  # Second part
```

### Be sure to set 
```
# email in which to receive job completion events
export JOB_MAIL=dude@cool.com

# Universal software path
export MUGQIC_INSTALL_HOME=/software/areas/genomics/phase2

```
# Hands on
We will follow the steps directly from the 2 steps.batch.[12].sh scripts

# Simpler way
We wrote pipelines that manage all that we did

This lists all the possbile processing steps

```
module load mugqic/mugqic_pipelines/2.1.0 mugqic/python/2.7.8
python $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.py -h
```

And this generates the script to run to actually do the processing

```
python $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.py -c rnaseq.course.ini -d design.tsv -r readsets.tsv -f -s 2,4-18 > steps.pbs.sh
```

