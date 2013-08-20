ART_454  README (last updated on 11/22/2011) Weichun Huang at whduke@gmail.com                              

DESCRIPTION:

	ART (art_454 version 1.6.8) is a simulation program to generate sequence read data of Roche 454
      	Pyrosequencing sequencers. ART generates reads according to the empirical read quality profile
       	and the calibrated error profile of uncall/overcall homopolymers from real 454 read data. ART 
	has been using for testing or benchmarking a variety of method or tools for next-generation 
	sequencing data analysis, including read alignment, de novo assembly, detection of SNP, CNV, or
       	other structure variation.

	art_454 can generate both single-end and paired-end of 454 sequencing platform.  Its outputs include
       	FASTQ read, ALN alignment, STAT read coverage, and optional SAM alignment files. ALN files can also
	be converted to UCSC BED files by using aln2bed.pl program included. 
       	
COMPILATION AND INSTALLATION:

	PREREQUISITES: 
			1) GNU g++ 4.0 or above (http://gcc.gnu.org/install) 
			2) GNU gsl library (http://www.gnu.org/s/gsl/) 

	./configure --prefix=$HOME
	make
	make install

EXAMPLES

	In the "examples" subdirectory, the shell script "run_test_examples.sh" gives two examples of using
       	art_454 for read simulation.  To test these two examples, just run the script "run_test_examples.sh"
USAGE

       SINGLE-END READ SIMULATION

          art_454 [ -p read_profile ] [-s] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <FOLD_COVERAGE>

          Example:
              art_454 -s seq_reference.fa ./outdir/dat_single_end 20

       PAIRED-END READ SIMULATION

          art_454 [ -p read_profile ] [-s] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <FOLD_COVERAGE> <MEAN_FRAG_LEN> <STD_DEV>

          Example:
             art_454 -s seq_reference.fa ./outdir/dat_paired_end 20 500 20

      COMMAND-LINE PARAMETERS

          INPUT_SEQ_FILE     -  the filename of DNA reference sequences in FASTA format
          OUTPUT_FILE_PREFIX -  the prefix or directory of output read data file (*.fq) and read alignment file (*.aln)
	  FOLD_COVERAGE      -  the fold of read coverage over the reference sequences 
	  MEAN_FRAG_LEN      -  the average DNA fragment size for paired-end read simulation 
	  STD_DEV            -  the standard deviation of the DNA fragment size for paired-end read simulation 

      OPTIONAL PARAMETERS 
      
         -s indicate to generate a SAM alignment file.   
         -p specify user's own read profile for simulation 
	 
	    the name of a read profile is the filename of directory containing read profile data files. a read profile
	    should consist of the four following files:  

	    1) qual_1st_profile   - empirical quality distribution of the 1st base of DNA homopolymers 
                                    These are homopolymer-length dependent distributions (one distribution for each length of homopolymer)
	    2) qual_mc_profile    - Markov Chain-based empirical quality distribution of the remaining bases of DNA homopolymers 
	    3) indel_erro_profile - undercall or overcall error profiles of DNA homopolymers 
	    4) length_dist        - the empirical distribution of 454 read length 

	    Please see the real example the under the folder "454_FLX_profile" for file formats of these four files. 

OUTPUT FILES

	*.fq   - FASTQ read data files. For paired-end read simulation, *1.fq contains data of the first reads, and *2.fq for the second reads.
	*.aln  - read alignment files. For paired-end read simulation, *1.aln has read alignments for the first reads and *2.aln for the second reads.
	*.sam  - SAM alignment file. 
	*.stat - read coverage files with the following three tab-delimited fields per lines: 
		 (1)reference position  
		 (2)number of reads starting at the position  
		 (3)number of reads covering the position  

	FASTQ FILE format 
		A FASTQ file contains both sequence bases and quality scores of sequencing reads and is in the following format:  
			@read_id 
			sequence_read 
			+ 
			base_quality_scores 
	
		A base quality score is coded by the ASCII code of a single character, where the quality score is equal to ASCII code of the
       		character minus 33.    

		Example: 
		@refid-4028550-1 
		caacgccactcagcaatgatcggtttattcacgat...
		+ 
		????????????7?????<??>??=&?<<?-<?0?...

	ALN format 
		An ALN file has a Header and main Body parts. The header part includes the command used to generate this file and reference
	       	sequence id and length. The header @CM tag for command line, and @SQ for reference sequence.  A header always starts with 
		"##ART" and ends with  "##Header End".

		HEADER EXAMPLE

		##ART_454
       		@CM     ../../bin/MacOS64/art_454 -s ./testSeq.fa ./single_454_test 10 
		@SQ     seq1    7207
	       	@SQ     seq2    3056
		##Header End

		The body part contains all alignments in the following format 

		>ref_seq_id	read_id	aln_start_pos	ref_seq_strand
	       	ref_seq_aligned
	       	read_seq_aligned 
	
		aln_start_pos is the alignment start position of reference sequence. aln_start_pos is always relative to the strand of reference
       		sequence. That is, aln_start_pos 10 in the plus (+) strand is different from aln_start_pos 10 in the minus (‐) stand.  
	
		ref_seq_aligned is the aligned region of reference sequence, which can be from plus strand or minus strand of the reference sequence.  
		read_seq_aligned is the aligned sequence read, which always in the same orientation of the same read in the corresponding fastq file. 

	SAM format

       		SAM is a standard format for next-gen sequencing read alignments. The details of the format and examples are available at the links below:	
		1) http://samtools.sourceforge.net/SAM1.pdf
	       	2) http://genome.sph.umich.edu/wiki/SAM		


	BED format

		See the format at UCSC http://genome.ucsc.edu/FAQ/FAQformat.html#format1

		NOTE: both ALN and BED format files use 0-based coordinate system while SAM format uses 1-based coordinate system.
