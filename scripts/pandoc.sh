#!/bin/sh

# IMPORTANT: always leave two empty spaces at the end of each md file,
# otherwise pandonc doesn't distinguish the sections properly
#
# human readable command:
# pandoc *.md
# -o ../pdf/practical.pdf
# --toc
# --variable title:"Cancer Genomics SNV, Exploratory practical - Cancer Genomics Workshop 2014"
# --variable date:"2014/07/02"
# --variable author:"Louis Letourneau"
# --variable links-as-notes
# --variable linkcolor:black
# --variable urlcolor:black
# --variable geometry:margin=3cm

# From the doc directory
pandoc *.md -o ../pdf/practical.pdf --toc --variable title:"Cancer Genomics SNV, Exploratory practical - Cancer Genomics Workshop 2014" --variable date:"2014/07/02" --variable author:"Louis Letourneau" --variable links-as-notes --variable linkcolor:black --variable urlcolor:black --variable geometry:margin=3cm


