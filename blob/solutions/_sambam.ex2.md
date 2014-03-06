Because all the reads from a same fragment always have the same names.
In this case since it's paired data, we should have 2 or more hits per read name.

I say or more because BWA mem can break reads into multiple pieces, so more than 2 hits can be found.
The use of the primary alignment flag or supplementary alignment flags can be used to identify these hits.
