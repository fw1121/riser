#!/bin/bash          

# location of fasta files 
str1=$1;
# read length
str2=$2;
# coverage
str3=$3;
# platform 
str4=$4;
# output directory
str5=$5;
# location of the executable
str6=$6;
# number of replicates
replicates=$7;

echo "###########################################   Simulate reads   #####################################################"
echo "Simulation will be run for the "${str4}" platform - with following parameters :";
echo "Read length: "${str2};
echo "Coverage: "${str3}"x";
echo "Output directory : "${str5};
echo "####################################################################################################################"

if [ "${str4}" == "illumina" ]; then
     p_type="art_illumina"
elif [ "${str4}" == "454" ]; then
     p_type="art_454"
elif [ "${str4}" == "solid" ]; then
     p_type="art_SOLiD"
fi

mut_fasta_location=${str1}
samples_full_path=($(ls -d -1 $mut_fasta_location/** | grep '.fa$'))
samples_names=($(ls $mut_fasta_location | grep '.fa$'))
resultsdir=${str5}/${str3}"xlen"${str2}"/fastq"
samples=${#samples_full_path[@]}

for ((i=0; i < $samples; i++ ));
#for i in 0;
do
  mkdir -p $resultsdir/${samples_names[$i]} 
  for (( j=0; j < $replicates; j++ ));
    do
    ${str6}/$p_type -i ${samples_full_path[$i]}  -o $resultsdir/${samples_names[$i]}/single_end.$j -l ${str2} -f ${str3} -sam
    # run calmd on the ART output 
    samtools calmd -S $resultsdir/${samples_names[$i]}/single_end.$j.sam ${samples_full_path[$i]} > $resultsdir/${samples_names[$i]}/single_end_md.$j.sam
    rm $resultsdir/${samples_names[$i]}/single_end.$j.sam
  done 
done
