#!/usr/bin/python
from __future__ import division
import re, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import itertools
from utils_make_txs import *
import random 
from random  import randint

random.seed(1234)

dir_in = sys.argv[1]
mutations = [int(x) for x in sys.argv[2].split(",")]

transcripts_in = dir_in + "/transcripts.txt"
fasta_in = dir_in + "/sequence.fa"

f = open(transcripts_in, 'r')
outfasta0 = open(dir_in + '/simulated_transcripts_0.fa', 'a')
outjunction = open(dir_in + '/simulated_junction_info.txt', 'a')

starts = {}
ends = {}
accession = str()
gene_name = {}
gene_id = {}
strand={}

seq_fasta = SeqIO.parse(open(fasta_in,"rU"), "fasta")

for record in seq_fasta:
        seq = record.seq 


print '---------------------------- Generate transcript sequences for the simulation --------------------------------------'
print '1. Reading the transcript file'


genecnt = 0

for line in f:
    genecnt+=1
    a = line.split('\t')    
    lst = a[9].split(',') 
    starts[a[0]] = [int(x) for x in lst]
    lst = a[10].rstrip('\n').split(',') 
    ends[a[0]] = [int(x) for x in lst]
    exonnmb = len(starts[a[0]])
    gene_name[a[0]] = a[1]
    gene_id[genecnt] = a[0]
    strand[a[0]] = a[3]
    if genecnt == 1:
        accession = a[2]
 
f.close()


print '2. Generate transcripts and junction info'


generate_transcripts(starts,ends,gene_name,strand,accession,seq,outjunction,outfasta0) 


print '3. Generate mutated sequences'
print '--------------------------------------------------------------------------------------------------------------------'



for mutrate in mutations:
    simulated_tx='/simulated_transcripts' + '_' + str(mutrate) + '.fa'
    outfasta = open(dir_in + simulated_tx, 'a')
    #print outfasta
    mseq=mutate_sequence(seq,mutrate)
    generate_transcripts_mut(starts,ends,gene_name,strand,accession,mseq,outfasta) 
    outfasta.close()


outfasta0.close()
outjunction.close()










