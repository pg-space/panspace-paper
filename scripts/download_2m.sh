# !/bin/bash

DIRSAVE=data/batches_bacteria_2m
BASEPATH=https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2

# 1. download _md5.txt with names of files
wget $BASEPATH/md5sum.txt

mkdir -p $DIRSAVE

# not complete-from-here
cut md5sum.txt -d" " -f3 | head -n 3 | \
while read f; \
do wget $BASEPATH/$f -O $DIRSAVE/$f;\
done
