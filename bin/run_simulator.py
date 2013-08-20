#!/usr/bin/python
import ConfigParser,os,sys
from subprocess import call
from glob import iglob
from shutil import move
from os.path import join

use_message = '''
Usage:
    python run_simulator.py
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def get_version():
    return "1.2.0"


def move_files(src_glob, dst_folder):
    for fname in iglob(src_glob):
        move(fname, join(dst_folder, os.path.basename(fname)))


def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1



def generate_transcripts_gbff(gbff_names,
                         gbff_locs,
                         dir_out,
                         gbff_mut_rates,
                         read_length,
                         fold_cov,
                         plaform_type,
                         sim_exec,
                         reps):
    #path=os.getcwd().rstrip('bin')
    path=os.getcwd()
    for id,name in enumerate(gbff_names):
        # first test if genbank file exists
        check_file=gbff_locs[id] + "/" + name + ".gbff" 
        err_file=call(["ls",check_file])
        if err_file != 0:
            print 'Error: GenBank in .ini file does not exist.' 
            sys.exit()
        gbff_dir_out_g = dir_out + name + "/genome"
        call(["mkdir", "-p", gbff_dir_out_g])
        exec1 = path + "/scripts/convert_genbank_to_transcript_table.pl" 
        call([exec1,check_file,gbff_dir_out_g,name])
        # make transcript-junction file
        exec1 = path + "/bin/make_transcripts.py"
        call([exec1,gbff_dir_out_g,gbff_mut_rates])
        fasta_dir=gbff_dir_out_g + "/fasta"
        call(["mkdir", "-p",fasta_dir])
        # clean up the directory
        move_files(gbff_dir_out_g + "/simulated_transcripts_*.fa", fasta_dir)
        # make dir for simulated reads 
        sim_out=dir_out + name + "/simulation/simulated/"
        call(["mkdir", "-p",sim_out])
        exec2 = path + "/scripts/run_simulation.sh" 
        # generate reads
        call([exec2,fasta_dir,read_length,fold_cov,plaform_type,sim_out,sim_exec,reps])



def generate_transcripts_custom(custom_names,
                         dir_out,
                         custom_mut_rates,
                         custom_fasta,
                         custom_txs,       
                         read_length,
                         fold_cov,
                         plaform_type,
                         sim_exec,
                         reps):
    #path=os.getcwd().rstrip('bin')
    path=os.getcwd()
    for id,name in enumerate(custom_names):
        # first test if genbank file exists
        err_file=call(["ls",custom_fasta])
        if err_file != 0:
            print 'Error: Fatsa in the .ini file does not exist.' 
            sys.exit()
        err_file=call(["ls",custom_txs])
        if err_file != 0:
            print 'Error: Transcripts info in the .ini file does not exist.' 
            sys.exit()            
        custom_dir_out_g = dir_out + name + "/genome" 
        call(["mkdir", "-p", custom_dir_out_g])
        lnk_name= custom_dir_out_g + '/' + 'transcripts.txt'
        call(["ln","-s",custom_txs,lnk_name])
        lnk_name= custom_dir_out_g + '/' + 'sequence.fa'
        call(["ln","-s",custom_fasta,lnk_name])
        exec1 = path + "/bin/make_transcripts_custom.py"
        call([exec1,custom_txs,custom_fasta,custom_mut_rates,custom_dir_out_g])
        fasta_dir = custom_dir_out_g + "/fasta"
        call(["mkdir", "-p",fasta_dir])
        # clean up the directory
        move_files(custom_dir_out_g + "/simulated_transcripts_*.fa", fasta_dir)
        # make dir for simulated reads 
        sim_out=dir_out + name + "/simulation/simulated/"
        call(["mkdir", "-p",sim_out])
        exec2 = path + "/scripts/run_simulation.sh" 
        # generate reads
        call([exec2,fasta_dir,read_length,fold_cov,plaform_type,sim_out,sim_exec,reps])


if len(sys.argv) > 1:
    if sys.argv[1] == "-v" or sys.argv[1] == "--version":
        print "RiSER v%s" % (get_version())
        exit(0)

# Read the .ini file

Config = ConfigParser.ConfigParser()
path=os.getcwd()
Config.read(path + "/config/config.ini")

gbff_names=ConfigSectionMap("genbank")['name'].split(",")
gbff_locs=ConfigSectionMap("genbank")['location'].split(",")
gbff_mut_rates=ConfigSectionMap("genbank")['mut_rates']
custom_names=ConfigSectionMap("custom")['name'].split(",")
custom_mut_rates=ConfigSectionMap("custom")['mut_rates']
custom_fasta=ConfigSectionMap("custom")['fasta']
custom_txs=ConfigSectionMap("custom")['transcripts']
project_dir=ConfigSectionMap("project_dir")['location']
project_name=ConfigSectionMap("project_dir")['name']
read_length=ConfigSectionMap("simulator")['read_length']
fold_cov=ConfigSectionMap("simulator")['fold_coverage']
plaform_type=ConfigSectionMap("simulator")['platform_type']
sim_exec=ConfigSectionMap("simulator")['executable']
reps=ConfigSectionMap("simulator")['replicates']


# Make the project directory

dir_out=project_dir + '/' + project_name + '/'
mkdir_err=call(["mkdir", "-p", dir_out])
if mkdir_err != 0:
        print 'Error: Could not create the project directory exiting.' 
        sys.exit()


# Generate transcripts and simulate reads 
print "RiSER v%s" % (get_version())
print "-------------------------------------------------------------------------------------------------------"
print "Simulation run."
print "-------------------------------------------------------------------------------------------------------"

if gbff_names[0] != 'None':
    generate_transcripts_gbff(gbff_names,gbff_locs,dir_out,gbff_mut_rates,read_length,fold_cov,plaform_type,sim_exec,reps)


if custom_names[0] != 'None':
    generate_transcripts_custom(custom_names,dir_out,custom_mut_rates,custom_fasta,custom_txs,read_length,fold_cov,plaform_type,sim_exec,reps)

