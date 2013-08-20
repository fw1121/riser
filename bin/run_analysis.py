#!/usr/bin/python
import ConfigParser,os,sys
from subprocess import call
from glob import iglob
from shutil import move
from os.path import join

use_message = '''
Usage:
    python run_analysis.py
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


def run_analysis(gbff_names,
                 aligners,
                 aligned_files,
                 project_dir,
                 cov_fold,
                 read_length):
    path=os.getcwd().rstrip('bin')
    for gbff_id,gbff_name in enumerate(gbff_names):
        for al_id,al_name in enumerate(aligners):
            path_sim_data=project_dir + '/' + gbff_name + '/' + 'simulation/simulated/' + cov_fold + 'x' + 'len' + read_length + '/fastq'
            junction_file=project_dir + '/' + gbff_name + '/genome/simulated_junction_info.txt'
            exec1 = path + "/scripts/run_analysis.sh"
            alignedx=aligned_files[gbff_id] + '/' + al_name
            call([exec1,path_sim_data,junction_file,alignedx,al_name,read_length])
        

if len(sys.argv) > 1:
    if sys.argv[1] == "-v" or sys.argv[1] == "--version":
        print "RiSER v%s" % (get_version())
        exit(0)


Config = ConfigParser.ConfigParser()
#Config.read("../config/config.ini")
path=os.getcwd()
Config.read(path + "/config/config.ini")


if ConfigSectionMap("genbank")['name'].split(",")[0] != 'None' and ConfigSectionMap("custom")['name'].split(",")[0] != 'None':
    gbff_names=ConfigSectionMap("genbank")['name'].split(",") + ConfigSectionMap("custom")['name'].split(",")
elif ConfigSectionMap("genbank")['name'].split(",")[0] == 'None':
    gbff_names=ConfigSectionMap("custom")['name'].split(",")
elif ConfigSectionMap("custom")['name'].split(",")[0] == 'None':
    gbff_names=ConfigSectionMap("genbank")['name'].split(",")

aligners=ConfigSectionMap("aligners")['name'].split(",")
aligned_files=ConfigSectionMap("aligners")['aligned_files'].split(",")

if len(gbff_names) != len(aligned_files):
    print 'Error Exiting: the # "names" in [genbank] differs from the # "aligned_files" in [aligners].' 
    exit(0)
 
project_dir=ConfigSectionMap("project_dir")['location'] + '/' + ConfigSectionMap("project_dir")['name']
read_length=ConfigSectionMap("simulator")['read_length']
cov_fold=ConfigSectionMap("simulator")['fold_coverage']

print "RiSER v%s" % (get_version()) 
print "-------------------------------------------------------------------------------------------------------"
print "Analysis run."
print "-------------------------------------------------------------------------------------------------------"


run_analysis(gbff_names,aligners,aligned_files,project_dir,cov_fold,read_length)

