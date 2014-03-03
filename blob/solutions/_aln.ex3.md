Mostly to save IO operations. BWA outputs a SAM, which is uncompressed text. This is a big file and takes up a lo of space.
Piping the output directly skips the SAM generation.

The problem though is that more RAM and processors are needed since we sort the output automatically.
One workaround is to just generate an unsorted BAM with the output using samtools.
One advantage here is that bwa mem outputs split reads in best alignment order.

Another option is to do this in 2 steps and pay the space cost.

