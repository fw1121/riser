#!/bin/bash
#454 test examples

art_454=../art_454
# 1) simulation of 454 single-end reads at 10X coverage using the built-in 454 FLX read  profile

$art_454 -s ./testSeq.fa ./single_454_test 10

#convert an aln file to a bed file
../aln2bed.pl single_454_test.bed single_454_test.aln

# 2) simulation of 454 paried-end reads at 5X coverage with the mean fragment size 500 
#    and standard deviation 20 using the built-in 454 FLX read  profile

$art_454 -s ./testSeq.fa ./paired_454_test 5 500 20

#convert both aln files to a bed file
../aln2bed.pl paired_454_test.bed paired_454_test1.aln paired_454_test2.aln



