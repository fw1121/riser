#!/usr/bin/python
from __future__ import division
import operator
import itertools
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import random 
from random  import randint
#from numpy import *


def generate_transcripts(starts,ends,gene_name,strand,accession,seq,outjunction,outfasta):
    """
    all transcripts will be generated as + strand  
    """
    for ids in gene_name:
        seqout = str()
        sumseqlength = 0
        junction = str()
        exoncnt=len(starts[ids])
        for eindx,start in enumerate(starts[ids]):
            end=ends[ids][eindx]
            seqout = seqout + seq[(start - 1):end]
            if exoncnt == 1:
               junction = '0'
               outfasta.write(">%s\n%s\n" % (
                       gene_name[ids] + '|' + ids + '|'+ accession + '|' + str(exoncnt),   
             seqout))
               string = gene_name[ids] + '|' + ids + '|'+ accession + '|' + str(exoncnt)
               outjunction.write(string + '\t' + strand[ids] + '\t' + junction + '\t' + str(start) + '\t' + str(end) + '\t' + '0' + '\n')
            elif exoncnt > 1 and eindx < (exoncnt  - 1):
               sumseqlength = sumseqlength + len(seq[(start - 1):end]) 
               junction = junction + str(sumseqlength) + ','
        if junction != '0' and exoncnt > 1:        
            junction = junction.rstrip('\,')
            outfasta.write(">%s\n%s\n" % (
                gene_name[ids] + '|' + ids + '|'+ accession + '|' + str(exoncnt),
             seqout))
            string = gene_name[ids] + '|' + ids + '|'+ accession + '|' + str(exoncnt) 
            starts_s=[str(x) for x in starts[ids]]
            ends_s=[str(x) for x in ends[ids]]
            outjunction.write(string + '\t' + strand[ids] + '\t' + junction + '\t' + ','.join(starts_s) + '\t' + ','.join(ends_s) + '\t' + '0' + '\n')  
 


def mutate_sequence(seq,mutrate):
    mseq=seq.tomutable() 
    bases=['A','G','C','T']
    nt_pos=sorted(random.sample(range(len(seq)),int(round(mutrate*len(seq)/100))))
    for nt in nt_pos:
        mseq[nt]=random.sample(list(set(bases) - set(mseq[nt])),int(1))[0]
    new_seq = mseq.toseq()
    return new_seq



def generate_transcripts_mut(starts,ends,gene_name,strand,accession,seq,outfasta):
    """
    all transcripts will be generated as + strand  
    """
    for ids in gene_name:
        seqout = str()
        sumseqlength = 0
        junction = str()
        exoncnt=len(starts[ids])
        for eindx,start in enumerate(starts[ids]):
            end=ends[ids][eindx]
            seqout = seqout + seq[(start - 1):end]
            if exoncnt == 1:
               outfasta.write(">%s\n%s\n" % (
                       gene_name[ids] + '|' + ids + '|'+ accession + '|' + str(exoncnt),   
             seqout))
        if exoncnt > 1:    
            outfasta.write(">%s\n%s\n" % (
                gene_name[ids] + '|' + ids + '|'+ accession + '|' + str(exoncnt),
             seqout))
 
