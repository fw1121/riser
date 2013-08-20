from utils_analysis import * 

def blast_to_cigar(hsp,readL):
    cigtuple=[]
    pseudo_cigar=[]
    mismatch_pos=[i for i, (left, right) in enumerate(zip(str(hsp.query),str(hsp.sbjct))) if left != right and left != '-' and right != '-']
    del_pos=[i for i, (left, right) in enumerate(zip(str(hsp.query),str(hsp.sbjct))) if left != right and left == '-' and right != '-']
    ins_pos=[i for i, (left, right) in enumerate(zip(str(hsp.query),str(hsp.sbjct))) if left != right and left != '-' and right == '-']
    soft_start=hsp.query_start - 1
    soft_end=readL - hsp.query_end
    for pos in range(0,len(str(hsp.query))):
        try:
            mismatch_pos.index(pos)       
            pseudo_cigar.append(8)   
        except ValueError:
            try:
                del_pos.index(pos)  
                pseudo_cigar.append(2)   
            except ValueError:
                try:
                    ins_pos.index(pos)  
                    pseudo_cigar.append(1)
                except ValueError:    
                    pseudo_cigar.append(7)             
    cigtuple = count(pseudo_cigar)
    if soft_start > 0:
       cigtuple = [(4,soft_start)] + cigtuple
    if soft_end > 0:
       cigtuple = cigtuple + [(4,soft_end)] 
    return cigtuple 



def base_accuracy_blast(cig_aln,cig_true,readSstart,readSend,coordinate_change):
    aln_stats_o = []
    aln_stats_no =[]
    refStart=int(coordinate_change[0])
    refEnd=int(coordinate_change[1])
    if not coordinate_change[3] or (coordinate_change[3] and coordinate_change[4]):
       range_t = range(refStart,refEnd+1)
       # case where you have clipped alignments at the beginning of the read  
       if cig_aln[0][0] == 4 or cig_aln[0][0] == 5:
           readSstart = readSstart - cig_aln[0][1] 
       range_aln = range(readSstart,readSend+1)
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
    elif coordinate_change[3] and not coordinate_change[4]:
        jonct_pos=sorted([coordinate_change[0],coordinate_change[1]] + list(coordinate_change[5])) 
        range_t=[]
        junction_lengths=[]
        for i in xrange(0,len(jonct_pos),2):
           range_t = range_t  + range(int(jonct_pos[i]),int(jonct_pos[i+1])+1)
           if i < (len(jonct_pos) - 2):
                junction_lengths = junction_lengths + [int(jonct_pos[i + 2]) - int(jonct_pos[i + 1]) - 1]
        # case where you have clipped alignments at the beginning of the read          
        if cig_aln[0][0] == 4 or cig_aln[0][0] == 5:
            readSstart = readSstart - cig_aln[0][1] 
        range_aln = range(readSstart,readSend+1)
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


