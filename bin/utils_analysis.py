from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam, re
from itertools import groupby
from numpy import *

def MD_tag(val):
    x = 0
    md = []
    for y in xrange(len(val)):
        if val[y].isalpha():
            if y == 0:
               offset = 0
	    elif y > 0 :   
              offset = int(val[x:y])
	    base = val[y]
	    md.append(offset)
	    md.append(base)
	    x = y + 1
    if x < len(val):
        md.append(int(val[x:]))
    return md



def find_index_all(L, value): 
    index = []
    i = -1
    try:
        while 1:
            i = L.index(value, i+1)
            index.append(i)
    except ValueError:
        pass
    return index
       

def count(data):
    tple = [(k, len(list(g))) for k, g in groupby(data)]
    return tple


def get_pseudo_seq(cig):
    pseudo_seq=[]
    for a in cig:
        pseudo_seq = pseudo_seq + [a[0]] * a[1]   
    return pseudo_seq


def parse_MD_tag(val):
    x = 0
    mdops = []
    del_split = val.split('^')
    qpseudo = [] 
    delition = ()
    if len(del_split) == 1:
      mdops = MD_tag(val)
    else: 
     delition = True
     nmb_del = len(del_split) - 1
     del_cnt = []
     val_ops = []
     if del_split[0] == '' or nmb_del > 1:
        del_cnt = range(1,nmb_del + 1)
        # all the other possibilities 
     elif del_split[0] != '' and nmb_del == 1:
       del_cnt = [1]
     for id,d in enumerate(del_split):
        if id in del_cnt:
           delstr = re.search(r'\A\D+', del_split[id]).group(0)    
           lendelstr = len(delstr)       
           delsub = str(lendelstr) + 'D'
           delsub = delsub + re.sub(r'^' + delstr, '', del_split[id])  
           val_ops.append(delsub)
        else:
           val_ops.append(d)
     for dels in val_ops:
         mdops = mdops + MD_tag(dels)
    if not delition:
        for mdraw in mdops:
          if str(mdraw).isalpha() and mdraw != '0':      
              qpseudo.append(int(8))
          elif str(mdraw).isdigit() and mdraw != '0': 
              for q in xrange(mdraw): 
                  qpseudo.append(int(7)) 
    elif delition:
          indxdel = find_index_all(mdops, 'D')   
          indxdelone = [x - 1 for x in indxdel]
          delcnt = 0
          cnt = 0 
          for mdraw in mdops:
            if str(mdraw).isalpha() and mdraw != '0' and mdraw != 'D': 
                qpseudo.append(int(8))
            elif mdraw == 'D':
                 for q in xrange(mdops[indxdelone[delcnt]]):
                      qpseudo.append(int(2))  
                 delcnt += 1  
            elif str(mdraw).isdigit() and str(mdraw) != '0' and not (cnt in indxdelone):
                 for q in xrange(mdraw):
                      qpseudo.append(int(7))
            cnt += 1
    return qpseudo 




def parse_cigar_MD_tags(cigar, MD):
        """
	Parse CIGAR and MD string and return pysam-suitable tuples.

	"""
        # parse the MD tag 
        if  isinstance(MD, basestring):       
            pseudoMD = parse_MD_tag(MD)        
        else:
            print 'Error: MD should be a string exiting.'
            #exit() 
        if isinstance(cigar, basestring):
           # first get the info from the cigar tag  
	   ops = dict( (i,j) for (j, i) in enumerate('MIDNSHP=X') )
	   cigtuple = []
	   n = ''
	   for c in cigar:
		if c.isdigit():
       			n += c
		elif c in ops:
			cigtuple.append( (ops[c], int(n)) )
			n = ''
        else: 
           cigtuple = cigar 
 	basesum = [0]
	cigtuplenew = []
        charsum = 0
        non_md = [1,3,4,5,6]
        id2 = 0
        for id,val in enumerate(cigtuple):
           if int(val[0]) in non_md:
              id2 += 1
           else:
              charsum = charsum + int(val[1]) 
	      basesum.append(charsum)       
           if int(val[0]) == 0 or int(val[0]) == 2:  
              subpseudomd  = pseudoMD[basesum[(id - id2)] : (basesum[(id - id2)] + int(val[1]))]
              mdtpl = count(subpseudomd)        
	      for tpl in mdtpl: 
                 cigtuplenew.append(tpl)
           else:
            cigtuplenew.append(val)
        return cigtuplenew   




def transfrom_MRNA_to_DNA_ref_frame(readTrueStart,readTrueEnd,exonE,exonS,exonj,readTrueLength,intronp):
        """
        transform coordinates from mRNA to DNA 
        ref frame
        """   
        junction_coord=[]
        cntjspan = 0
        cntjunction = 0
        junction_type = str()
        junction = False 
        exon_len = 0
        exon_len_array = []
        nmbexons = range(len(exonS))
        cnt_exon = 0
        adjoin=False
        junctions_hit=[] 
        for jindx in exonj:
            cntjunction += 1
            if int(readTrueStart) <=  int(jindx) and int(readTrueEnd) > int(jindx):
               cntjspan += 1
               junction = True
               junctions_hit.append(jindx)
               if (int(exonS[exonj.index(jindx) + 1]) - 1) == (int(exonE[exonj.index(jindx)])):  
                   adjoin=True  
               if cntjunction == intronp and intronp != 0 :
                  junction_type = 'EI'   
               else:
                  junction_type = 'EE'
        for exonindx in nmbexons: 
               exon_len =  exon_len + int(exonE[exonindx]) - int(exonS[exonindx]) + 1
               exon_len_array.append(exon_len)    
	       if exon_len_array[exonindx] > readTrueStart and cnt_exon  == 0 and not junction:
                  DNArefStart =  int(exonS[exonindx]) + int(readTrueStart) - 1
                  DNArefEnd = int(DNArefStart) + readTrueLength - 1
                  break
	       elif exon_len_array[exonindx] >= readTrueStart and cnt_exon  == 0 and junction:
                  DNArefStart =  int(exonS[exonindx]) + int(readTrueStart) - 1
		  DNArefEnd = int(exonS[(cnt_exon + 1)]) + (readTrueLength - (int(exonE[exonindx]) - int(DNArefStart))) - 2
                  junction_coord=(int(exonE[exonindx]),int(exonS[(cnt_exon + 1)]))
		  break
	       elif(exon_len_array[exonindx] > readTrueStart and cnt_exon  >= 1 and not junction):  
                  DNArefStart =  (int(exonS[exonindx]) + int(readTrueStart) - exon_len_array[(cnt_exon - 1)]) - 1
                  DNArefEnd = int(DNArefStart) + readTrueLength - 1
                  break
	       elif(exon_len_array[exonindx] >= readTrueStart and cnt_exon  >= 1 and junction): 
                  DNArefStart =  (int(exonS[exonindx]) + int(readTrueStart) - exon_len_array[(cnt_exon - 1)]) - 1
                  if cnt_exon < len(exonS) - 1:
                     DNArefEnd = int(exonS[(cnt_exon + 1)]) + (readTrueLength - (int(exonE[exonindx]) - int(DNArefStart))) - 2
                     junction_coord=(int(exonE[exonindx]),int(exonS[(cnt_exon + 1)])) 
                  elif cnt_exon ==  len(exonS) - 1:
                     # in case read crossing junction/last exon  
                     DNArefEnd = int(exonS[(cnt_exon)]) + (readTrueLength - (int(exonS[(cnt_exon)]) - int(DNArefStart))) - 1
                     junction_coord=(int(exonE[cnt_exon-1]),int(exonS[(cnt_exon)])) 
		  break
               cnt_exon += 1
        # case of reads spanning multiple junctions 
        if cntjspan > 1:
            exon_len_s = 0
            exon_len_array_s = []
            for exonindx_s in nmbexons: 
               exon_len_s =  exon_len_s + int(exonE[exonindx_s]) - int(exonS[exonindx_s]) + 1
               exon_len_array_s.append(exon_len_s)
            j_start=exonj.index(junctions_hit[0])     
            DNArefEnd = int(exonS[j_start + cntjspan]) + (int(readTrueStart) + int(readTrueLength) - exon_len_array_s[j_start + cntjspan - 1]) - 1 - 1
            junction_coord=[]
            for i,j in enumerate(exonE[(j_start):(j_start + cntjspan)]):
                junction_coord=junction_coord + [int(j),int(exonS[j_start + 1:(j_start + cntjspan + 1)][i])]
        return DNArefStart,DNArefEnd,junction_type,junction,adjoin,junction_coord,cntjspan





def base_accuracy(cig_aln,cig_true,readSstart,readSend,coordinate_change):
    aln_stats_o = []
    aln_stats_no =[]
    refStart=int(coordinate_change[0])
    refEnd=int(coordinate_change[1])
    if not coordinate_change[3] or (coordinate_change[3] and coordinate_change[4] and coordinate_change[6] == 1):
       range_t = range(refStart,refEnd+1)
       # case where you have clipped alignments at the beginning of the read  
       if cig_aln[0][0] == 4 or cig_aln[0][0] ==5:
           readSstart = readSstart - cig_aln[0][1] 
       readSend=readSstart + cig_aln[:,1].sum() - cig_aln[:,1][cig_aln[:,0] == 1].sum()
       range_aln = range(readSstart,readSend)
       cig_aln_p_o=[]
       cig_true_p_o=[]
       cig_aln_p_o_ins=[]
       cig_true_p_o_ins=[]
       cig_aln_p=get_pseudo_seq(cig_aln)
       cig_true_p=get_pseudo_seq(cig_true)
       # position of gaps along the read 
       gap_pos_aln=find_index_all(cig_aln_p, 1)
       gap_pos_true=find_index_all(cig_true_p, 1)
       gap_pos_overlap=list(set(gap_pos_true) & set(gap_pos_aln))
       for x in sorted(list(set(range_t) & set(range_aln))):
           if len(set([range_aln.index(x)]).intersection(set(gap_pos_aln))) == 0: 
               cig_aln_p_o.append(cig_aln_p[range_aln.index(x)])
           if len(set([range_t.index(x)]).intersection(set(gap_pos_true))) == 0:    
               cig_true_p_o.append(cig_true_p[range_t.index(x)])
       base_accu=len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left == right])
       aln_stats_o.append((1,len(gap_pos_aln),len(gap_pos_true),len(gap_pos_overlap)))
       aln_stats_o.append((2,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 2])))
       aln_stats_o.append((3,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 3])))
       aln_stats_o.append((4,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 4])))
       aln_stats_o.append((5,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 5])))
       aln_stats_o.append((6,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 6])))
       aln_stats_o.append((7,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 7])))
       aln_stats_o.append((8,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 8])))
       aln_stats_no.append((1,len(gap_pos_aln),len(gap_pos_true),len(gap_pos_overlap)))
       aln_stats_no.append((2,cig_aln_p.count(2) - cig_aln_p_o.count(2)))
       aln_stats_no.append((3,cig_aln_p.count(3) - cig_aln_p_o.count(3)))
       aln_stats_no.append((4,cig_aln_p.count(4) - cig_aln_p_o.count(4)))
       aln_stats_no.append((5,cig_aln_p.count(5) - cig_aln_p_o.count(5)))                   
       aln_stats_no.append((6,cig_aln_p.count(6) - cig_aln_p_o.count(6)))  
       aln_stats_no.append((7,cig_aln_p.count(7) - cig_aln_p_o.count(7)))
       aln_stats_no.append((8,cig_aln_p.count(8) - cig_aln_p_o.count(8)))
    elif (coordinate_change[3] and not coordinate_change[4] and coordinate_change[6] == 1) or (coordinate_change[3] and coordinate_change[6] > 1):
        jonct_pos=sorted([coordinate_change[0],coordinate_change[1]] + list(coordinate_change[5])) 
        range_t=[]
        junction_lengths=[]
        for i in xrange(0,len(jonct_pos),2):
           range_t = range_t  + range(int(jonct_pos[i]),int(jonct_pos[i+1])+1)
           if i < (len(jonct_pos) - 2):
                junction_lengths = junction_lengths + [int(jonct_pos[i + 2]) - int(jonct_pos[i + 1]) - 1]
        # case where you have clipped alignments at the beginning of the read          
        if cig_aln[0][0] == 4 or cig_aln[0][0] ==5:
            readSstart = readSstart - cig_aln[0][1]
        readSend=readSstart + cig_aln[:,1].sum() - cig_aln[:,1][cig_aln[:,0] == 1].sum()
        range_aln = range(readSstart,readSend)
        cig_aln_p_o=[]
        cig_true_p_o=[]
        cig_aln_p_o_ins=[]
        cig_true_p_o_ins=[]
        cig_aln_p=get_pseudo_seq(cig_aln)
        cig_true_p=get_pseudo_seq(cig_true)
        # position of gaps along the read 
        gap_pos_aln=find_index_all(cig_aln_p, 1)
        gap_pos_true=find_index_all(cig_true_p, 1)
        gap_pos_overlap=list(set(gap_pos_true) & set(gap_pos_aln))
        for x in sorted(list(set(range_t) & set(range_aln))):
            if len(set([range_aln.index(x)]).intersection(set(gap_pos_aln))) == 0:
                cig_aln_p_o.append(cig_aln_p[range_aln.index(x)])
            if len(set([range_t.index(x)]).intersection(set(gap_pos_true))) == 0:
                cig_true_p_o.append(cig_true_p[range_t.index(x)])
        base_accu=len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left == right])
        aln_stats_o.append((1,len(gap_pos_aln),len(gap_pos_true),len(gap_pos_overlap)))
        aln_stats_o.append((2,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 2])))
        aln_stats_o.append((3,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 3])))
        aln_stats_o.append((4,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 4])))
        aln_stats_o.append((5,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 5])))
        aln_stats_o.append((6,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 6])))
        aln_stats_o.append((7,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 7])))
        aln_stats_o.append((8,len([i for i, (left, right) in enumerate(zip(cig_true_p_o,cig_aln_p_o)) if left != right and  right == 8])))
        aln_stats_no.append((1,len(gap_pos_aln),len(gap_pos_true),len(gap_pos_overlap)))
        aln_stats_no.append((2,cig_aln_p.count(2) - cig_aln_p_o.count(2)))
        aln_stats_no.append((3,cig_aln_p.count(3) - cig_aln_p_o.count(3)))
        aln_stats_no.append((4,cig_aln_p.count(4) - cig_aln_p_o.count(4)))
        aln_stats_no.append((5,cig_aln_p.count(5) - cig_aln_p_o.count(5)))                   
        aln_stats_no.append((6,cig_aln_p.count(6) - cig_aln_p_o.count(6)))  
        aln_stats_no.append((7,cig_aln_p.count(7) - cig_aln_p_o.count(7)))
        aln_stats_no.append((8,cig_aln_p.count(8) - cig_aln_p_o.count(8)))
    return base_accu,aln_stats_o,aln_stats_no


