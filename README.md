riser
=====
**RiSER** 

Downloading and using **RiSER** is free, if you use **RiSER** or its code in your work please acknowledge **RiSER** by referring to its GitHub homepage https://github.com/oicr-ibc/riser

This is important for us since obtaining grants is one significant way to fund
planning and implementation for our projects. Also if you find **RiSER** useful
in your research feel free to let us know.

RiSER is brought to you by:
 * Vincent Ferretti
 * Ivan Borozan 
 * Stuart Watt


RiSER was originally developed by:
 * Ivan Borozan


Getting Started
---------------

Minimum Requirements
-----------------------
Tested on UBUNTU-12.04

__R (2.14.1):__
```bash
$ apt-get install r-base
```

Python (2.7.3) - the program assumes that Python is in `/usr/bin/python`
perl (5.14.2)
samtools (0.1.18)

Following Perl module needs to be installed:

__bioperl:__
```bash
$ apt-get install bioperl
```

Following Python modules need to be installed in the order shown below:

If you do not have pip installed, install it as shown below:
```bash
$ sudo apt-get install python-pip python-dev  
```

__numpy(1.6.2):__
```bash
$ sudo pip install numpy
```

__BioPython(1.6):__
```bash
$ sudo pip install biopython 
```

__rpy2(2.3.1):__
```bash
$ sudo pip install rpy2
```

__setuptools(1.0):__
```bash
$ sudo easy_install -U distribute
```

__Cython(0.17.4):__
```bash
$ sudo pip install cython
```

__pysam(0.6):__
```bash
$ sudo pip install pysam
```

Installation
------------
This version of RiSER has has been tested under Linux (Ubuntu 12.04).

To install:

Option 1:
    
```bash
$ sudo apt-get install git
$ git clone https://github.com/oicr-ibc/riser.git 
$ cd riser
```

Option 2:

```bash
$ wget https://github.com/oicr-ibc/riser/archive/master.zip
$ unzip master.zip
$ mv riser-master riser
$ cd riser
```


Done! No installation is required, all Python scripts are in `$RISER_DIR/bin` and should be compatible with your system

Usage: 	 
------  

Assuming you are working in the RiSER directory. 	

1. Simulate data for a particular set of genomes. 

    The default `config.ini` file is in the `config` directory - to run RiSER you need to modify the initial `.ini` file. However make sure to first run RiSER with the configuration file `config/config_simulation.example` provided as an example on how to run the simulation:

    ```bash
    $ cp config/config.ini config/config.ini.save
    $ cp config/config_simulation.example config/config.ini
    ```

    edit the `config.ini` file.

    run the simulation script:

    ```bash    
    $ python bin/run_simulator.py
    ```

    (i) Note, simulated results will be output to the directory specified in the `config.ini` file (see `config_simulation.example`). Also make sure to check if you are running a 32-bit or 64-bit machine (see `config_simulation.example` under `[aligners]`)

    (ii) Note, in the GenBank flat file, the GenBank 'FEATURES' entries 'gene' and 'CDS' if both present, need to have the /db_xref="GeneID:XXXXX" associated with each. 


2. Compare the aligner's output to the truth file (from simulated data):

    Make sure to first run RiSER with the configuration file `config/config_analysis.example` provided as an example on how to run the analysis:

    ```bash
    $ cp config/config.ini config/config.ini.save
    $ cp config/config_analysis.example config/config.ini
    ```

    edit the `config.ini`	
	
    run the analysis script:

    ```bash
    $ python bin/run_analysis.py
    ```

    The summary of results will be output to the aligner's directory specified by the user in the `config.ini` file (see `config_analysis.example`)

    In the example provided, results for the NC_001357.1 genome and the BFAST aligner will be output to:

    `examples/aligners/NC_001357.1/BFAST/simulated_transcripts_0.fa/Rdata_multi/`

    `examples/aligners/NC_001357.1/BFAST/simulated_transcripts_10.fa/Rdata_multi/`
   	
    Note that the results are output as R data files, to view them launch R and load results as shown below: 

    ```bash
    $ cd examples/aligners/NC_001357.1/BFAST/simulated_transcripts_0.fa/Rdata_multi/
    ```

    in R type:

    ```r
    > # To load the data
    > load("aligner_stats.gzip")
    > 
    > # To run the analysis statistics
    > aligner_stats
    ```
   
**File format for user specified transcript files**:	 

In the case a transcript file is specified by the user (see also `examples/genomes/NC_001357.1_transcripts.txt`) each row in the file should designate a single transcript and columns (tab delimited) should be set as in the order shown below:

transcript_id (e.g. GI number or any other unique id) \t transcript_name \t genome_id (e.g. GenBank Accession) \t strand \t transcript_START \t transcript_END \t transcript_START \t transcript_END \t numb_exons \t exons_START(the START positions of each exon needs to be separated by commas) \t exons_END(the END positions of each exon needs to be separated by commas (and in the same order as the START positions))   	  	             						       

Datasets
--------

More datasets are available on the wiki at: https://github.com/oicr-ibc/riser/wiki/Datasets.


License and Copyright
---------------------
Licensed under the GNU General Public License, Version 3.0. See LICENSE for more details.

Copyright 2013 The Ontario Institute for Cancer Research.

Acknowledgement
---------------
This project is supported by the Ontario Institute for Cancer Research
(OICR) through funding provided by the government of Ontario, Canada.
