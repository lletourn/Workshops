#---PREPARING THE REFERENCE -----

cd /lb/project/mugqic/projects/cshl201411/
for i in G_*;do mkdir -p ${i}/browser/draft ; mkdir -p ${i}/browser/miseq ; mkdir -p ${i}/browser/final;done

#Assembly = multi fasta file of un-ono contigs CONTIG_FILE.FA
#To select the contigs >1kb:
# - put the fasta file on one line
for i in G_*;do 
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/draft
  ln -s ../../pacbio_assembly/${i}/30X/merSize14/report/consensus.fasta.gz ./
  zcat consensus.fasta.gz | perl /lb/project/mugqic/projects/cshl201411/clean_newline_multifasta.pl /dev/stdin > ${i}.draft.fasta
done

for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/miseq
  cat Scaffolds.fasta | perl /lb/project/mugqic/projects/cshl201411/clean_newline_multifasta.pl /dev/stdin > ${i}.miseq.fasta
done

for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/final
  perl /lb/project/mugqic/projects/cshl201411/clean_newline_multifasta.pl polish1_out/data/consensus.fasta > ${i}.final.fasta
done


#Rename contigs part 1-X
cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do 
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/draft
  # order contigs from biggest to smallest
  perl /lb/project/compgen/jwassers/Perl/sort_fasta_by_length.pl ${i}.draft.fasta | awk 'BEGIN {IDX=1} {if($0 ~ /^>/){print ">part" IDX;IDX+=1;} else {print $0}}' > ${i}.draft.sorted.fasta

  #(#Extra step for bacterial genomes (circular): we want to cut the 1st contig in the middle (call the pieces contigN_start and contigN_end), move the 1st piece at the end of the fasta file > circ.CONTIG_FILE.FA)
 
  # To insert 10kb spacers in between contigs
  perl  /lb/project/mugqic/projects/cshl201411/add_spacers_in_ref.pl  ${i}.draft.sorted.fasta  10000  | perl /lb/project/mugqic/projects/cshl201411/clean_newline_multifasta.pl /dev/stdin > ${i}.draft.ref.fa
done

# MiSeq
cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/miseq;
  perl /lb/project/compgen/jwassers/Perl/sort_fasta_by_length.pl  ${i}.miseq.fasta | awk 'BEGIN {IDX=1} {if($0 ~ /^>/){print ">part"  IDX;IDX+=1;} else {print $0}}' > ${i}.miseq.sorted.fasta;

  perl  /lb/project/mugqic/projects/cshl201411/add_spacers_in_ref.pl   ${i}.miseq.sorted.fasta  30000  | perl   /lb/project/mugqic/projects/cshl201411/clean_newline_multifasta.pl   /dev/stdin > ${i}.miseq.ref.fa;
done

# Final, circularized
cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/final;
  perl /lb/project/compgen/jwassers/Perl/sort_fasta_by_length.pl  ${i}.final.fasta | awk 'BEGIN {IDX=1} {if($0 ~ /^>/){print ">part"  IDX;IDX+=1;} else {print $0}}' > ${i}.final.sorted.fasta ; 
  perl  /lb/project/mugqic/projects/cshl201411/add_spacers_in_ref.pl   ${i}.final.sorted.fasta 30000  | perl   /lb/project/mugqic/projects/cshl201411/clean_newline_multifasta.pl   /dev/stdin > ${i}.final.ref.fa;
done

#---PREPARING THE CONTIG TRACK-----
 
#Blat the contigs against the ref, select best hit and do pslToBed
module load mugqic/ucsc/20140212
cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  echo $i;
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/draft;
  blat ${i}.draft.ref.fa ${i}.draft.sorted.fasta -stepSize=5 -tileSize=11 -minScore=0 -minIdentity=0 -noHead contigs.psl;
  awk  '{OFS="\t"} {print $1+$3-$2-$5-$7, $0}' contigs.psl |  sort -S1G +10 -11 +0 -1nr | awk '!($11 in l){print;l[$11]=1}' | cut -f 2-22 > bh.contigs.psl;
  pslToBed bh.contigs.psl bh.contigs.bed;
done

cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  echo $i;
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/miseq;
  blat ${i}.miseq.ref.fa ${i}.miseq.sorted.fasta -stepSize=5 -tileSize=11 -minScore=0 -minIdentity=0 -noHead contigs.psl;
  awk  '{OFS="\t"} {print $1+$3-$2-$5-$7, $0}' contigs.psl |  sort -S1G +10 -11 +0 -1nr | awk '!($11 in l){print;l[$11]=1}' | cut -f 2-22 > bh.contigs.psl;
  pslToBed bh.contigs.psl bh.contigs.bed;
done

cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  echo $i;
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/final;
  blat ${i}.final.ref.fa ${i}.final.sorted.fasta -stepSize=5 -tileSize=11 -minScore=0 -minIdentity=0 -noHead contigs.psl;
  awk  '{OFS="\t"} {print $1+$3-$2-$5-$7, $0}' contigs.psl |  sort -S1G +10 -11 +0 -1nr | awk '!($11 in l){print;l[$11]=1}' | cut -f 2-22 > bh.contigs.psl;
  pslToBed bh.contigs.psl bh.contigs.bed;
done
#---TO CHANGE THE TRACK DEFINITIONS and UPDATE THE BROWSER -----
 
#The location of your track configuration file :
/data/share/ucsc/track25112010/Leishmania/BUILD/trackDb.ra
 
#Whenever that file gets modified, to "refresh" the browser :
cd /data/share/ucsc/track25112010
make alpha DBS=BUILD

#---PACBIO DATASET---------*
module load mugqic/blast/2.2.29+

cd /lb/project/mugqic/projects/cshl201411/
for i in G*;do
  echo $i;
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/;
  ln -s ../pacbio_assembly/${i}/filtering/data/filtered_subreads.fasta ./;

  #Filter 3kb+ PacBio reads
  perl /lb/project/mugqic/projects/cshl201411/clean_newline_multifasta.pl filtered_subreads.fasta > filtered_subreads.oneLine.fasta;
  perl -p -e "s/>(\S+)\n/>\1 /g;" filtered_subreads.oneLine.fasta | awk '{if(length($2)>=3000){print}}' | perl -p -e "s/>(\S+) />\1\n/g;"  > 3kb.PB.fa;
  perl -p -e "s/>(\S+)\n/>\1 /g;" filtered_subreads.oneLine.fasta | awk '{if(length($2)>=7000){print}}' | perl -p -e "s/>(\S+) />\1\n/g;"  > 7kb.PB.fa;

  #Create blastable database
  cd draft;
  makeblastdb -in ${i}.draft.ref.fa -dbtype nucl;

  #Blast PacBio reads
  blastn -task dc-megablast -query ../3kb.PB.fa -strand both -db ${i}.draft.ref.fa  -outfmt 6 -out 3kb.PB.out &
  blastn -task dc-megablast -query ../7kb.PB.fa -strand both -db ${i}.draft.ref.fa  -outfmt 6 -out 7kb.PB.out &
done

# Miseq
cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/miseq;
  makeblastdb -in ${i}.miseq.ref.fa -dbtype nucl;
  blastn -task dc-megablast -query ../3kb.PB.fa -strand both -db ${i}.miseq.ref.fa  -outfmt 6 -out 3kb.PB.out &
  blastn -task dc-megablast -query ../7kb.PB.fa -strand both -db ${i}.miseq.ref.fa  -outfmt 6 -out 7kb.PB.out &
done

# Final
cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/final;
  makeblastdb -in ${i}.final.ref.fa -dbtype nucl;
  blastn -task dc-megablast -query ../3kb.PB.fa -strand both -db ${i}.final.ref.fa  -outfmt 6 -out 3kb.PB.out &
  blastn -task dc-megablast -query ../7kb.PB.fa -strand both -db ${i}.final.ref.fa  -outfmt 6 -out 7kb.PB.out &
done

cd /lb/project/mugqic/projects/cshl201411/
module load mugqic/exonerate;
for i in G_[46]*;do
  echo $i
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/draft/;
  #Keep hits for criterias 80% id + 500bp span min
  awk '{print $0,$8-$7}' 3kb.PB.out | awk '{if ($3>=80 && $13>=500) print $0}' > 80.500.3kb.PB.out;
  awk '{print $0,$8-$7}' 7kb.PB.out | awk '{if ($3>=80 && $13>=500) print $0}' > 80.500.7kb.PB.out;

  #Select best hit (bh)
  sort -k 13nr -k 1 80.500.3kb.PB.out | awk '{if(l[$1]<1){print;l[$1]+=1}}' > bh.80.500.3kb.PB.out;
  sort -k 13nr -k 1 80.500.7kb.PB.out | awk '{if(l[$1]<1){print;l[$1]+=1}}' > bh.80.500.7kb.PB.out;

  #Make PacBio bed track
  fastalength ${i}.draft.ref.fa;
  SZ=`fastalength ${i}.draft.ref.fa | cut -f1 -d' '`
  perl /lb/project/compgen/jwassers/Perl/generate_tracks_for_pacbio_reads.pl bh.80.500.3kb.PB.out ${SZ} >  bh.80.500.3kb.PB.bed;
  perl /lb/project/compgen/jwassers/Perl/generate_tracks_for_pacbio_reads.pl bh.80.500.7kb.PB.out ${SZ} >  bh.80.500.7kb.PB.bed;
done

# Miseq
cd /lb/project/mugqic/projects/cshl201411/
module load mugqic/exonerate;
for i in G_*;do
  echo $i
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/miseq/;
  awk '{print $0,$8-$7}' 3kb.PB.out | awk '{if ($3>=80 && $13>=500) print $0}' > 80.500.3kb.PB.out;
  awk '{print $0,$8-$7}' 7kb.PB.out | awk '{if ($3>=80 && $13>=500) print $0}' > 80.500.7kb.PB.out;

  sort -k 13nr -k 1 80.500.3kb.PB.out | awk '{if(l[$1]<1){print;l[$1]+=1}}' > bh.80.500.3kb.PB.out;
  sort -k 13nr -k 1 80.500.7kb.PB.out | awk '{if(l[$1]<1){print;l[$1]+=1}}' > bh.80.500.7kb.PB.out;
  fastalength ${i}.miseq.ref.fa;
  SZ=`fastalength ${i}.miseq.ref.fa | cut -f1 -d' '`
  perl /lb/project/compgen/jwassers/Perl/generate_tracks_for_pacbio_reads.pl bh.80.500.3kb.PB.out ${SZ} >  bh.80.500.3kb.PB.bed;
  perl /lb/project/compgen/jwassers/Perl/generate_tracks_for_pacbio_reads.pl bh.80.500.7kb.PB.out ${SZ} >  bh.80.500.7kb.PB.bed;
done

# Final
module load mugqic/exonerate;
cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  echo $i
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/final;
  awk '{print $0,$8-$7}' 3kb.PB.out | awk '{if ($3>=80 && $13>=500) print $0}' > 80.500.3kb.PB.out;
  awk '{print $0,$8-$7}' 7kb.PB.out | awk '{if ($3>=80 && $13>=500) print $0}' > 80.500.7kb.PB.out;

  sort -k 13nr -k 1 80.500.3kb.PB.out | awk '{if(l[$1]<1){print;l[$1]+=1}}' > bh.80.500.3kb.PB.out;
  sort -k 13nr -k 1 80.500.7kb.PB.out | awk '{if(l[$1]<1){print;l[$1]+=1}}' > bh.80.500.7kb.PB.out;
  fastalength ${i}.final.ref.fa;
  SZ=`fastalength ${i}.final.ref.fa | cut -f1 -d' '`
  perl /lb/project/compgen/jwassers/Perl/generate_tracks_for_pacbio_reads.pl bh.80.500.3kb.PB.out ${SZ} >  bh.80.500.3kb.PB.bed;
  perl /lb/project/compgen/jwassers/Perl/generate_tracks_for_pacbio_reads.pl bh.80.500.7kb.PB.out ${SZ} >  bh.80.500.7kb.PB.bed;
done

*#---ILLUMINA DATASET----------*
# Assemble with Ray
cd G_089/browser/miseq/ray/k31
# try k21, k31, k41
vim ray.pbs # adjust kmer et al.
qsub ray.pbs

# after assembly, launch flanking k to find best asm.

# convert the fastq in fasta
cd /lb/project/mugqic/projects/cshl201411/
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/GROUP1-089_S5_L001_R1_001.fastq.gz miseq_cshlws_10282014/M02986/GROUP1-089_S5_L001_R1_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_089/browser/MISEQ_1.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/GROUP1-089_S5_L001_R2_001.fastq.gz miseq_cshlws_10282014/M02986/GROUP1-089_S5_L001_R2_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_089/browser/MISEQ_2.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/GROUP2-6920_S6_L001_R1_001.fastq.gz miseq_cshlws_10282014/M02986/GROUP2-6920_S6_L001_R1_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_6920/browser/MISEQ_1.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/GROUP2-6920_S6_L001_R2_001.fastq.gz miseq_cshlws_10282014/M02986/GROUP2-6920_S6_L001_R2_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_6920/browser/MISEQ_2.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/GROUP3-6919_S7_L001_R1_001.fastq.gz miseq_cshlws_10282014/M02986/GROUP3-6919_S7_L001_R1_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_6914/browser/MISEQ_1.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/GROUP3-6919_S7_L001_R2_001.fastq.gz miseq_cshlws_10282014/M02986/GROUP3-6919_S7_L001_R2_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_6914/browser/MISEQ_2.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/GROUP4-4681_S8_L001_R1_001.fastq.gz miseq_cshlws_10282014/M02986/GROUP4-4681_S8_L001_R1_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_4681/browser/MISEQ_1.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/GROUP4-4681_S8_L001_R2_001.fastq.gz miseq_cshlws_10282014/M02986/GROUP4-4681_S8_L001_R2_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_4681/browser/MISEQ_2.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/VM-690_S1_L001_R1_001.fastq.gz miseq_cshlws_10282014/M02986/VM-690_S1_L001_R1_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_690/browser/MISEQ_1.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/VM-690_S1_L001_R2_001.fastq.gz miseq_cshlws_10282014/M02986/VM-690_S1_L001_R2_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_690/browser/MISEQ_2.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/VM-698_S2_L001_R1_001.fastq.gz miseq_cshlws_10282014/M02986/VM-698_S2_L001_R1_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_698/browser/MISEQ_1.fasta &
zgrep -A1 "^@M0[02]" miseq_cshlws_10282014/M00557/VM-698_S2_L001_R2_001.fastq.gz miseq_cshlws_10282014/M02986/VM-698_S2_L001_R2_001.fastq.gz | grep -v "\-\-"|sed 's/.*gz://g'|perl -p -e "s/\@M0/\>M0/;" > G_698/browser/MISEQ_2.fasta &
time wait


cd /lb/project/mugqic/projects/cshl201411/
module load mugqic/ucsc/20140212 mugqic/samtools/0.1.19-gpfs mugqic/bedtools/2.17.0
for i in G_[46]*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/;
  cat MISEQ_1.fasta MISEQ_2.fasta > MISEQ.fasta;

  # keep only the 1st 150bp if necessary
  #   cut -b 1-150 MISEQ.fasta > MiSeq.150.fa

  #Blat reads
  cd draft;
  blat ${i}.draft.ref.fa ../MISEQ.fasta -stepSize=5 -tileSize=11 -minScore=0 -minIdentity=0 -noHead MiSeq.psl &
done;
time wait

# MiSeq
module load mugqic/ucsc/20140212 mugqic/samtools/0.1.19-gpfs mugqic/bedtools/2.17.0
cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/miseq;
  blat ${i}.miseq.ref.fa ../MISEQ.fasta -stepSize=5 -tileSize=11 -minScore=0 -minIdentity=0 -noHead MiSeq.psl &
done;
time wait

# Final
module load mugqic/ucsc/20140212 mugqic/samtools/0.1.19-gpfs mugqic/bedtools/2.17.0
cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/final;
  blat ${i}.final.ref.fa ../MISEQ.fasta -stepSize=5 -tileSize=11 -minScore=0 -minIdentity=0 -noHead MiSeq.psl &
done;
time wait

cd /lb/project/mugqic/projects/cshl201411/
for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/draft;
  #Select all perfect hits (ph)
  awk  '{if($1==$11 && $2==0 && $3==0 && $4==0 && $5==0 && $6==0 && $7==0 && $8==0){print $0}}' MiSeq.psl > allph.MiSeq.psl
  awk  '{if($1 >= ($11*0.8) && $2==0 && $3==0 && $4==0 && $5==0 && $6==0 && $7==0 && $8==0){print $0}}' MiSeq.psl > allGood.MiSeq.psl

  #Make bed track
  pslToBed allph.MiSeq.psl allph.MiSeq.bed
  pslToBed allGood.MiSeq.psl allGood.MiSeq.bed

  #Make 0_coverage_in_MiSeq_reads track
  samtools faidx ${i}.draft.ref.fa
  genomeCoverageBed -d -i allph.MiSeq.bed -g ${i}.draft.ref.fa.fai > CVG.TBL

  #  select positions where coverage=0
  awk '{if($3==0){print $0;}}' CVG.TBL > 0CVG

  #  make the bed track
  perl /lb/project/mugqic/projects/cshl201411/generate_track_of_covered_genome.pl 0CVG > NO_CVG_IN_MISEQ.bed
done

# MiSeq
cd /lb/project/mugqic/projects/cshl201411/
module load mugqic/ucsc/20140212 mugqic/samtools/0.1.19-gpfs mugqic/bedtools/2.17.0
for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/miseq;
  awk  '{if($1==$11 && $2==0 && $3==0 && $4==0 && $5==0 && $6==0 && $7==0 && $8==0){print $0}}' MiSeq.psl > allph.MiSeq.psl
  awk  '{if($1 >= ($11*0.8) && $2==0 && $3==0 && $4==0 && $5==0 && $6==0 && $7==0 && $8==0){print $0}}' MiSeq.psl > allGood.MiSeq.psl
  pslToBed allph.MiSeq.psl allph.MiSeq.bed
  pslToBed allGood.MiSeq.psl allGood.MiSeq.bed
  samtools faidx ${i}.miseq.ref.fa
  genomeCoverageBed -d -i allph.MiSeq.bed -g ${i}.miseq.ref.fa.fai > CVG.TBL
  awk '{if($3==0){print $0;}}' CVG.TBL > 0CVG
  perl /lb/project/mugqic/projects/cshl201411/generate_track_of_covered_genome.pl 0CVG > NO_CVG_IN_MISEQ.bed
done

# Final
module load mugqic/ucsc/20140212 mugqic/samtools/0.1.19-gpfs mugqic/bedtools/2.17.0
cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/final;
  awk  '{if($1==$11 && $2==0 && $3==0 && $4==0 && $5==0 && $6==0 && $7==0 && $8==0){print $0}}' MiSeq.psl > allph.MiSeq.psl
  awk  '{if($1 >= ($11*0.8) && $2==0 && $3==0 && $4==0 && $5==0 && $6==0 && $7==0 && $8==0){print $0}}' MiSeq.psl > allGood.MiSeq.psl
  pslToBed allph.MiSeq.psl allph.MiSeq.bed
  pslToBed allGood.MiSeq.psl allGood.MiSeq.bed
  samtools faidx ${i}.final.ref.fa
  genomeCoverageBed -d -i allph.MiSeq.bed -g ${i}.final.ref.fa.fai > CVG.TBL
  awk '{if($3==0){print $0;}}' CVG.TBL > 0CVG
  perl /lb/project/mugqic/projects/cshl201411/generate_track_of_covered_genome.pl 0CVG > NO_CVG_IN_MISEQ.bed
done

cd /lb/project/mugqic/projects/cshl201411/ ; for i in G_*;do
  cd /lb/project/mugqic/projects/cshl201411/${i}/browser/final;
  zgrep -v "^#" motif_out/data/motifs.gff.gz | perl -n -e 'my @val=split(/\t/); $val[8] =~ /.*IPDRatio=([-.0-9]+)/ ; if($val[6] eq "+") {print "chr1\t".($val[3]-1)."\t".$val[4]."\t$1\n"}' >ipdratio_F.bed
  zgrep -v "^#" motif_out/data/motifs.gff.gz | perl -n -e 'my @val=split(/\t/); $val[8] =~ /.*IPDRatio=([-.0-9]+)/ ; if($val[6] eq "-") {print "chr1\t".($val[3]-1)."\t".$val[4]."\t$1\n"}' >ipdratio_R.bed
done


#---PREPARING THE BROWSER------
On the vervet server
 
export CLIENT_NAME=CSHL_2014
/data/share/bin/cgb/cgb create_browser
/data/share/bin/cgb/cgb add_clade       CSHL_2014 CSHL_2014 1
/data/share/bin/cgb/cgb add_genome G_089 CSHL_2014 1
/data/share/bin/cgb/cgb add_genome G_690 CSHL_2014 2
/data/share/bin/cgb/cgb add_genome G_698 CSHL_2014 3
/data/share/bin/cgb/cgb add_genome G_4681 CSHL_2014 4
/data/share/bin/cgb/cgb add_genome G_6914 CSHL_2014 5
/data/share/bin/cgb/cgb add_genome G_6920 CSHL_2014 6

cd /data/share/lletourn/cshl201411/;
for i in G_[46]*;do
  /data/share/bin/cgb/cgb add_build ${i}_draft "Pacbio Draft build" ${i} part1 ${i} source 1234

  #Load the reference into the browser
  /data/share/bin/cgb/cgb add_fasta        ${i}_draft part1 ${i}/browser/draft/${i}.draft.ref.fa
  /data/share/bin/cgb/cgb add_defaultdb    ${i} ${i}_draft
done

cd /data/share/lletourn/cshl201411/;for i in G_*;do
  /data/share/bin/cgb/cgb add_build ${i}_miseq "MiSeq Build" ${i} part1 ${i} source 1234

  /data/share/bin/cgb/cgb add_fasta        ${i}_miseq part1 ${i}/browser/miseq/${i}.miseq.ref.fa
done

cd /data/share/lletourn/cshl201411/;for i in G_*;do
  /data/share/bin/cgb/cgb add_build ${i}_final "Final build" ${i} part1 ${i} source 1234
  /data/share/bin/cgb/cgb add_fasta        ${i}_final part1 ${i}/browser/final/${i}.final.ref.fa
done

#Turn on Blat for that reference :
# check which ports are available :
ps -afe | grep gfServer

# pick 2 ports that aren't used yet (random 5 digit numbers i.e 58100 58101)
/data/share/bin/cgb/cgb add_blat  G_089_draft  50100  50101
/data/share/bin/cgb/cgb add_blat  G_690_draft  50102  50103
/data/share/bin/cgb/cgb add_blat  G_698_draft  50104  50105
/data/share/bin/cgb/cgb add_blat  G_4681_draft  50106  50107
/data/share/bin/cgb/cgb add_blat  G_6914_draft  50108  50109
/data/share/bin/cgb/cgb add_blat  G_6920_draft  50110  50111

/data/share/bin/cgb/cgb add_blat  G_089_miseq  50200  50201
/data/share/bin/cgb/cgb add_blat  G_690_miseq  50202  50203
/data/share/bin/cgb/cgb add_blat  G_698_miseq  50204  50205
/data/share/bin/cgb/cgb add_blat  G_4681_miseq  50206  50207
/data/share/bin/cgb/cgb add_blat  G_6914_miseq  50208  50209
/data/share/bin/cgb/cgb add_blat  G_6920_miseq  50210  50211


/data/share/bin/cgb/cgb add_blat  G_089_final  50300  50301
/data/share/bin/cgb/cgb add_blat  G_690_final  50302  50303
/data/share/bin/cgb/cgb add_blat  G_698_final  50304  50305
/data/share/bin/cgb/cgb add_blat  G_4681_final  50306  50307
/data/share/bin/cgb/cgb add_blat  G_6914_final  50308  50309
/data/share/bin/cgb/cgb add_blat  G_6920_miseq  50310  50311

cd /data/share/lletourn/cshl201411/;for i in G_*;do
  cd /data/share/lletourn/cshl201411/${i}/browser/draft/;
  #Load the contig track in the browser
  hgLoadBed ${i}_draft contigs bh.contigs.bed;

  #Load the Pacbio track
  hgLoadBed ${i}_draft Pacbio3kb bh.80.500.3kb.PB.bed
  hgLoadBed ${i}_draft Pacbio7kb bh.80.500.7kb.PB.bed

  #Load the track Illumina
  hgLoadBed ${i}_draft MiSeq_all_perfect_hits allph.MiSeq.bed
  hgLoadBed ${i}_draft MiSeq_all_good_hits  allGood.MiSeq.psl

  #Load the NoCVG track
  hgLoadBed ${i}_draft no_cvg_in_MiSeq NO_CVG_IN_MISEQ.bed
done

cd /data/share/lletourn/cshl201411/;for i in G_*;do
  cd /data/share/lletourn/cshl201411/${i}/browser/miseq;
  hgLoadBed ${i}_miseq contigs bh.contigs.bed;
  hgLoadBed ${i}_miseq Pacbio3kb bh.80.500.3kb.PB.bed
  hgLoadBed ${i}_miseq Pacbio7kb bh.80.500.7kb.PB.bed
  hgLoadBed ${i}_miseq MiSeq_all_perfect_hits allph.MiSeq.bed
  hgLoadBed ${i}_miseq MiSeq_all_good_hits allGood.MiSeq.bed
  hgLoadBed ${i}_miseq no_cvg_in_MiSeq NO_CVG_IN_MISEQ.bed
done

cd /data/share/lletourn/cshl201411/;for i in G_*;do
  cd /data/share/lletourn/cshl201411/${i}/browser/final ; 
  hgLoadBed ${i}_final contigs bh.contigs.bed;
  hgLoadBed ${i}_final Pacbio3kb bh.80.500.3kb.PB.bed
  hgLoadBed ${i}_final Pacbio7kb bh.80.500.7kb.PB.bed
  hgLoadBed ${i}_final MiSeq_all_perfect_hits allph.MiSeq.bed
  hgLoadBed ${i}_final MiSeq_all_good_hits allGood.MiSeq.bed
  hgLoadBed ${i}_final no_cvg_in_MiSeq NO_CVG_IN_MISEQ.bed

  hgLoadBed ${i}_final ipdratio_F ipdratio_F.bed
  hgLoadBed ${i}_final ipdratio_R ipdratio_R.bed
done

# Update trackDb.ra
cp /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_draft/trackDb.ra /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_miseq/trackDb.ra ; cp /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_draft/trackDb.ra /data/share/ucsc/trackDb20140214/CSHL_2014/G_698_miseq/trackDb.ra ; cp /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_draft/trackDb.ra /data/share/ucsc/trackDb20140214/CSHL_2014/G_4681_miseq/trackDb.ra ; cp /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_draft/trackDb.ra /data/share/ucsc/trackDb20140214/CSHL_2014/G_6914_miseq/trackDb.ra ; cp /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_draft/trackDb.ra /data/share/ucsc/trackDb20140214/CSHL_2014/G_6920_miseq/trackDb.ra

cd /data/share/ucsc/trackDb20140214
make alpha DBS=G_089_miseq ; make alpha DBS=G_690_miseq ; make alpha DBS=G_698_miseq ; make alpha DBS=G_4681_miseq ; make alpha DBS=G_6914_miseq ; make alpha DBS=G_6920_miseq

cp /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_final/trackDb.ra /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_final/trackDb.ra ; cp /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_final/trackDb.ra /data/share/ucsc/trackDb20140214/CSHL_2014/G_698_final/trackDb.ra ; cp /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_final/trackDb.ra /data/share/ucsc/trackDb20140214/CSHL_2014/G_4681_final/trackDb.ra ; cp /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_final/trackDb.ra /data/share/ucsc/trackDb20140214/CSHL_2014/G_6914_final/trackDb.ra ; cp /data/share/ucsc/trackDb20140214/CSHL_2014/G_089_final/trackDb.ra /data/share/ucsc/trackDb20140214/CSHL_2014/G_6920_final/trackDb.ra

cd /data/share/ucsc/trackDb20140214
make alpha DBS=G_089_final ; make alpha DBS=G_690_final ; make alpha DBS=G_698_final ; make alpha DBS=G_4681_final ; make alpha DBS=G_6914_final ; make alpha DBS=G_6920_final
