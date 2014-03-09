#!/bin/sh

# IMPORTANT: always leave two empty spaces at the end of each md file,
# otherwise pandonc doesn't distinguish the sections properly
#
# human readable command:
# pandoc *.md
# -o ../pdf/practical.pdf
# --toc
# --variable title:"DNA-Seq processing - Kyoto Course on Bioinformatics 2014"
# --variable date:"2014/03/10"
# --variable author:"Louis Letourneau"
# --variable links-as-notes
# --variable linkcolor:black
# --variable urlcolor:black
# --variable geometry:margin=3cm

pandoc *.md -o ../pdf/practical.pdf --toc --variable title:"DNA-Seq processing - Kyoto Course on Bioinformatics 2014" --variable date:"2014/03/10" --variable author:"Louis Letourneau" --variable links-as-notes --variable linkcolor:black --variable urlcolor:black --variable geometry:margin=3cm


