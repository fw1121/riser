#!/bin/bash 

# dir for simulated data
dir_data_simu=$1

# junction_file
junction_file=$2

# dir for aligned data
dir_aligned_data=$3

# aligner
aligner=$4

# read length
readL=$5

dirs_simu=($(ls $dir_data_simu ))
dirs_nmb_simu=${#dirs_simu[@]}
dirs_aln=($(ls -d $dir_aligned_data/* | grep '.fa'))
dirs_nmb_aln=${#dirs_aln[@]}

if [ $aligner = "BLAST" ]; then
 for ((i=0; i < $dirs_nmb_simu; i++));
 #for i in 0;
 do
     echo $dir_data_simu/${dirs_simu[$i]}/
     #echo ${dirs_aln[$i]}/sam/bam
     echo $junction_file
     mkdir ${dirs_aln[$i]}/Rdata_multi
     python ./bin/analysis_BLAST.py $dir_data_simu/${dirs_simu[$i]}/ ${dirs_aln[$i]}/xml_out/ $junction_file ${dirs_aln[$i]} $readL
 done
else
  for ((i=0; i < $dirs_nmb_simu; i++));
  #for i in 0;
 do
     echo $dir_data_simu/${dirs_simu[$i]}/
     echo ${dirs_aln[$i]}/sam/bam/
     echo $junction_file
     mkdir ${dirs_aln[$i]}/Rdata_multi	
     python ./bin/analysis_aligners.py $dir_data_simu/${dirs_simu[$i]}/ ${dirs_aln[$i]}/sam/bam/ $junction_file ${dirs_aln[$i]}
 done
fi
