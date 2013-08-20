from __future__ import division
from Bio.Blast import NCBIXML
from utils_parse_BLAST import * 
from utils_analysis import * 
import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from numpy import *
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2ri
import os, os.path
import sys


directory_simulated_data=sys.argv[1]
simulated_files=[]

for name in os.listdir(directory_simulated_data):
    files = re.search('.sam',name)
    if files is not None:
        simulated_files.append(directory_simulated_data + name)


simulated_files=sorted(simulated_files)


directory_aligner_data=sys.argv[2]
aligner_files=[]


for name in os.listdir(directory_aligner_data):
    files = re.search('.xml',name)
    if files is not None:
        aligner_files.append(directory_aligner_data + name)

aligner_files=sorted(aligner_files)


if len(simulated_files) != len(aligner_files):
   print 'Error: files should be of same size.'
   exit()


# read the junction info file 
junction_file = open(sys.argv[3], 'r')

t_strand={}
exon_j_pos={}
exonStart={}
exonEnd={}
intron_pos={}

for line in junction_file:
    jdata = line.split('\t')    
    t_strand[jdata[0]] = jdata[1]
    exon_j_pos[jdata[0]] = jdata[2].split(',')   
    exonStart[jdata[0]] = jdata[3].split(',')   
    exonEnd[jdata[0]] = jdata[4].split(',')   
    intron_pos[jdata[0]] = jdata[5].rstrip('\n').split(',')

junction_file.close()

tot_reads=[]
tot_reads_jctn=[]
mapped_aln_tot_non_jctn=[]
mapped_aln_tot_jctn=[]
base_accu_non_jctn=[]
base_accu_jctn=[]
stats_aln_mean_non_jctn=[]
stats_aln_mean_jctn=[]
stats_aln_mean_non_jctn2=[]
stats_aln_mean_jctn2=[]

for index00, file00 in enumerate(simulated_files):
 print 'reading file:', index00, file00

 samfile_true = pysam.Samfile(file00, 'r')

 readTrueStart = []
 readTrueEnd = []
 readTrueId = []
 readTrueRef = []
 readTrueStrand = []
 readTrueSeq = []
 readTrueCig = []
 readTrueLength = []
 readMDTrue = []
 scoresTrue=[]
 avgQualTrue=[]
 low_base_quality=[]

 for align in samfile_true.fetch():
    scoresTrue = [ord(m)-33 for m in align.qqual]
    scores_a = array(scoresTrue)
    # python base 0 need to add 1 to the start position
    readTrueStart.append(int(align.pos)  + 1)
    align_length = align.rlen
    readTrueLength.append(align.rlen)
    read_end = int(align.pos) + align_length 
    readTrueEnd.append(read_end)
    readTrueId.append(align.qname)
    readTrueRef.append(samfile_true.getrname(align.tid))
    strand = -1 if align.is_reverse else 1
    readTrueStrand.append(strand)
    readTrueSeq.append(align.query)
    readTrueCig.append(align.cigar)
    try:
      readMDTrue.append(align.opt('MD'))
    except KeyError:
      readMDTrue.append(None)    

 # read the aligner file 
 result_handle = open(aligner_files[index00])
 blast_records = NCBIXML.parse(result_handle)
 print 'Read aligner file:', aligner_files[index00]
 readstart = []
 readend = []
 readid = []
 readRef = []
 readstrand = []
 readseq = []
 readcig = []
 readAS=[]
 readqstart=[]
 readqend=[]

 for blast_record in blast_records:
    for alignment in blast_record.alignments:
       for hsp in alignment.hsps:
         readstrand.append(hsp.frame[1])
         if hsp.frame[1] == -1:
            readstart.append(hsp.sbjct_end)    
            readend.append(hsp.sbjct_start)
         elif hsp.frame[1] == 1:
            readstart.append(hsp.sbjct_start)   
            readend.append(hsp.sbjct_end) 
         readqstart.append(hsp.query_start)
         readqend.append(hsp.query_end) 
         readid.append(str(blast_record.query)) 
         readRef.append(str(alignment.title))
         readcig.append(blast_to_cigar(hsp,int(sys.argv[5])))
         readAS.append(hsp.bits)

 read_junction=[]
 junction_type=[]
 junction_adjoint=[]
 stats_o=[]
 stats_no=[]
 base_accu=[]
 strand_accuracy=[]
 junction_type_lstats=[]
 junction_adjoint_lstats=[]
 mapped_read=[]
 for read in readTrueId:
    # in case of multiple alignments 
    # readid.index(read) will get the first alignment which is the most sig one
    try:
        indx_aln = find_index_all(readid,read)
        if len(indx_aln) == 0:  
            mapped_read.append(False)          
            indx_aln=''  
        else:
            mapped_read.append(True)          
    except ValueError:
        indx_aln = ''
        mapped_read.append(False)          
    indx = readTrueId.index(read)
    coordinate_change = transfrom_MRNA_to_DNA_ref_frame(readTrueStart[indx],readTrueEnd[indx],exonEnd[readTrueRef[indx]],exonStart[readTrueRef[indx]],exon_j_pos[readTrueRef[indx]],readTrueLength[indx],int(intron_pos[readTrueRef[indx]][0])) 
    if len(indx_aln) == 1:
        indx_aln = indx_aln[0]
    elif len(indx_aln) > 1:
        max_AS_score = max([readAS[int(x)] for x in indx_aln])   
        AS_scores_a = array([readAS[int(x)] for x in indx_aln])
        indx_aln_a = array(indx_aln)
        for indx_AS in indx_aln:
            if (readAS[int(indx_AS)] == max_AS_score) and (readstart[indx_AS] - int(coordinate_change[0])) == 0:
                # take the read with max AS and best position  
                indx_aln = int(indx_AS)
                break 
            else:
                # take the first read with max AS 
                indx_aln = indx_aln_a[AS_scores_a == max_AS_score][0] 
    if indx_aln != '':
        read_junction.append(coordinate_change[3])
        junction_adjoint.append(coordinate_change[4])
        if  readstrand[indx_aln] == readTrueStrand[indx]:
            strand_accuracy.append(True)
        else:
            strand_accuracy.append(False)     
        if coordinate_change[3]:
            junction_type.append(coordinate_change[2])
        else:
            junction_type.append('')
        cig_aln = array(readcig[indx_aln]) 
        if readstrand[indx_aln] == readTrueStrand[indx]:
            if readstrand[indx_aln] == -1:
                cig_aln = cig_aln[::-1] 
            cig_true = array(readTrueCig[indx])
            stats=base_accuracy_blast(cig_aln,cig_true,readstart[indx_aln],readend[indx_aln],coordinate_change)
            base_accu.append(stats[0])
            stats_o.append(stats[1])
            stats_no.append(stats[2])
        else:
            base_accu.append(0)
            stats_o.append([(1,0,0,0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (7, 0), (8, 0)])
            stats_no.append([(1,0,0,0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (7, 0), (8, 0)])
    else:   
        read_junction.append(coordinate_change[3])
        if coordinate_change[3]:
            junction_type.append(coordinate_change[2])
            junction_adjoint.append(coordinate_change[4])
        else:
            junction_type.append('') 
            junction_adjoint.append(True) 


 print 'Runnning calc on stats.'

 # Generate stats:
 # total number of simulated reads
 tot_nmb_reads=len(readTrueId)
 # total number of simulated junctions 
 jt=array(junction_type)
 j_adj=array(junction_adjoint)
 j_adj_ee=j_adj[jt == 'EE']
 jt_ee=jt[jt == 'EE']
 tot_nmb_reads_jctn=len(jt_ee[j_adj_ee == False])
 tot_reads.append(tot_nmb_reads)
 tot_reads_jctn.append(tot_nmb_reads_jctn)

 # total numb of mapped reads across non-junctions
 strand_accuracy_a=array(strand_accuracy)
 mapped_read_a=array(mapped_read)
 jt_mapped=jt[mapped_read_a == True][strand_accuracy_a == True]
 mapped_aln_tot_non_jctn.append(len(jt_mapped[jt_mapped == ''])/(tot_nmb_reads - tot_nmb_reads_jctn)*100)

 
 if len(mapped_read_a[mapped_read_a == True]) == 0 or len(base_accu) == 0:
        mapped_aln_tot_jctn.append(0)
        base_accu_non_jctn.append(0)
        base_accu_jctn.append(0)
        stats_aln_mean_jctn.append([0,0,0,0,0,0,0,0])
        stats_aln_mean_non_jctn.append([0,0,0,0,0,0,0,0])
        stats_aln_mean_non_jctn2.append([0,0,0,0,0,0,0,0])
        stats_aln_mean_jctn2.append([0,0,0,0,0,0,0,0])
        continue 

 # total number of mapped reads across junctions
 j_adj_mapped=j_adj[mapped_read_a == True][strand_accuracy_a == True]
 j_adj_mapped_ee=j_adj_mapped[jt_mapped == 'EE']
 jt_mapped_ee=jt_mapped[jt_mapped == 'EE']
 mapped_aln_tot_jctn.append(len(jt_mapped_ee[j_adj_mapped_ee == False])/tot_nmb_reads_jctn*100)

 # base accuracy across non-junctions
 base_accu_a=array(base_accu)
 base_accu_non_jctn.append(base_accu_a[strand_accuracy_a == True][jt_mapped == ''].mean())

 if mapped_aln_tot_jctn[index00] == 0:
    base_accu_jctn.append(0)
 else:
    base_accu_jctn.append(base_accu_a[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False].mean()) 


 del base_accu_a,jt,j_adj,mapped_read_a

 # calculate stats on alignments
 # inside of alignment 
 stats_aln=[]
 stats_o_a = array(stats_o)

 for x in stats_o_a:
    stats_aln.append([x[0][1] - x[0][3], x[1][1], x[2][1], x[3][1], x[4][1], x[5][1], x[6][1], x[7][1]])

 stats_aln_a=array(stats_aln)

 # stats on alignment across non_junctions
 stats_aln_mean_non_jctn.append([stats_aln_a[strand_accuracy_a == True][jt_mapped == ''][:,0].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == ''][:,1].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == ''][:,2].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == ''][:,3].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == ''][:,4].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == ''][:,5].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == ''][:,6].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == ''][:,7].mean()])


 if mapped_aln_tot_jctn[index00] == 0:
     stats_aln_mean_jctn.append([0,0,0,0,0,0,0,0])
 else:
     stats_aln_mean_jctn.append([stats_aln_a[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,0].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,1].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,2].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,3].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,4].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,5].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,6].mean(),stats_aln_a[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,7].mean()])
 del stats_aln_a



 # outside of alignment 
 stats_aln2=[]
 stats_no_a = array(stats_no)

 for x in stats_no_a:
    stats_aln2.append([x[0][1] - x[0][3], x[1][1], x[2][1], x[3][1], x[4][1], x[5][1], x[6][1], x[7][1]])

 stats_aln_a2=array(stats_aln2)

 # stats on alignment across non_junctions
 stats_aln_mean_non_jctn2.append([stats_aln_a2[strand_accuracy_a == True][jt_mapped == ''][:,0].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == ''][:,1].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == ''][:,2].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == ''][:,3].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == ''][:,4].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == ''][:,5].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == ''][:,6].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == ''][:,7].mean()])

 # stats on alignment across junctions
 if mapped_aln_tot_jctn[index00] == 0:
     stats_aln_mean_jctn2.append([0,0,0,0,0,0,0,0])
 else:
     stats_aln_mean_jctn2.append([stats_aln_a2[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,0].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,1].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,2].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,3].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,4].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,5].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,6].mean(),stats_aln_a2[strand_accuracy_a == True][jt_mapped == 'EE'][j_adj_mapped_ee == False][:,7].mean()])

 del stats_aln_a2

# get stats out as R data file 

r_out_file=sys.argv[4] + '/Rdata_multi/aligner_stats.gzip'
stats1 = r['cbind'](tot_reads,tot_reads_jctn,mapped_aln_tot_non_jctn,mapped_aln_tot_jctn,base_accu_non_jctn,base_accu_jctn)
stats1.colnames=ro.StrVector(['Total reads(non-jctn)','Total reads(jctn)','Aligned reads(non-jctn)','Aligned reads(jctn)','Base Accu(non-jctn)','Base Accu(jctn)']) 
r.assign("aligner_stats", stats1)
r("save(aligner_stats," + "file=" "'" + r_out_file + "'"  + ", compress=TRUE)")


r_out_file=sys.argv[4] + '/Rdata_multi/stats_o.gzip'
stats2=r['cbind'](stats_aln_mean_non_jctn,stats_aln_mean_jctn)
stats2.colnames=ro.StrVector(['Non-jctn','Jctn'])
robj = numpy2ri(stats2)
r.assign("stats_o", robj)
r("save(stats_o," + "file=" "'" + r_out_file + "'"  + ", compress=TRUE)")


r_out_file=sys.argv[4] + '/Rdata_multi/stats_no.gzip'
stats3=r['cbind'](stats_aln_mean_non_jctn2,stats_aln_mean_jctn2)
stats3.colnames=ro.StrVector(['Non-jctn','Jctn'])
robj = numpy2ri(stats3)
r.assign("stats_no", robj)
r("save(stats_no," + "file=" "'" + r_out_file + "'"  + ", compress=TRUE)")
























