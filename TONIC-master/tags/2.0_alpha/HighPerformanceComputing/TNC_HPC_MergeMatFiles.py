#! /usr/bin/python                        
#
# Copyright (C) 2013 by Howard Hughes Medical Institute.
#
# Purpose: submit Matlab feature extraction job(s) to cluster

# ------------------------- imports -------------------------

import os, glob, gc
import shutil, commands
import sys, re, optparse

tonic_home = "/groups/zlatic/home/denisovg/work/Josh_Dudman/TONIC_v2" # default location for Sun Grid processing
tonic_data = "/groups/zlatic/home/denisovg/work/Josh_Dudman" 
if "TONIC_HOME" in os.environ.keys():
    tonic_home = os.environ['TONIC_HOME']
if "TONIC_DATA" in os.environ.keys():
    tonic_data = os.environ['TONIC_DATA']

# ----------------------------------------------------------------------

def file_merging_command_line_parser(parser):
    parser.add_option("-D", "--debug",     dest="debug", help="debugging mode; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-s", "--submission_command", dest="submission_command", help="source, qsub or qsub_debug", metavar="submission_command", default="qsub")
    parser.add_option("-v", "--verbose",   action="store_true", dest="verbose",               help="increase the verbosity level of output", default=False)
    parser.add_option("-c", "--compile",   action="store_true", dest="compile", help="compile Matlab code", metavar="compile", default=False)
    parser.add_option("-e", "--executable_name", dest="executable_name", help="executable", metavar="executable_name", default="TNC_HPC_MergeMatFiles")
    parser.add_option("-d", "--output_folder", dest="output_folder", help="folder where to place executable", metavar="output_folder", default=tonic_home + "/..")
    return parser

# -----------------------------------------------------------------------

def compile_code(options):
    input_path      = os.path.join(tonic_home,"HighPerformanceComputing", \
                                   "TNC_HPC_MergeMatFiles.m")
    command = "mcc -m " + input_path + " -o " + options.executable_name + \
              " -d " + options.output_folder
    if options.verbose:
        command += " -v "
        print "command= ", command, "\n"
    os.system(command)
    os.system("chmod +x " + os.path.join(options.output_folder, \
                                         options.executable_name))

# ----------------------------------------------------------------------

def get_maximal_indices(outfolderpath, prefix, my_type):
    max_shank = -1
    max_seg = -1
    for my_file in os.listdir(outfolderpath):
        if  re.search(prefix,  my_file) and \
            re.search("shank", my_file) and \
            re.search("seg",   my_file) and \
            re.search("_" + my_type, my_file):
            print "my_file=", my_file
            shank_ind =  int(my_file.split("shank")[1].split("_")[0])
            seg_ind   =  int(my_file.split("seg")[1].split("_")[0])
            if  max_shank < shank_ind:
                max_shank = shank_ind 
            if  max_seg   < seg_ind:   
                max_seg   = seg_ind    
    return (max_shank, max_seg)

# ----------------------------------------------------------------------

def create_shell_script(input_file_path, my_type, shank, options):
    
    outfolderpath = tonic_data
    # command1 = merge files
    # command2 = remove files
#   scr.write("export MCR_CACHE_ROOT=/scratch/$USER/mcr_cache_root.$JOB_ID.$SGE_TASK_ID\n")
#   scr.write("mkdir -p $MCR_CACHE_ROOT\n")
    command1 = os.path.join(tonic_home, options.output_folder, \
                            options.executable_name) + " "
    command2 = "rm -f "
    tokens = input_file_path.split("/")
    input_file = tokens[len(tokens)-1]
    prefix = input_file.split(".")[0]

    max_shank, max_seg = get_maximal_indices(outfolderpath, prefix, my_type)
    output_file_path = os.path.join(outfolderpath, prefix + "_" + my_type + ".mat")
    for k in range(1, max_seg+1):
        for j in range(1, max_shank+1):
            input_mat_file_path = \
                os.path.join(outfolderpath, prefix + "_shank" + str(j) + \
                             "_seg" + str(k) + "_" + my_type + ".mat")
            command1 += input_mat_file_path + " "
            command2 += input_mat_file_path + " "
    command1 += output_file_path

    shell_script_path = os.path.join(outfolderpath, "Shell_script_merge." +\
                                     prefix + "." + str(my_type) + ".sh")
    scr = open(shell_script_path, 'wt')
    scr.write("#!/usr/bash\n")
    scr.write(command1 + "\n")
    if not options.debug:
        scr.write(command2 + "\n")
    command1 += "\n"
    print "shell script command1=", command1, \
          "             command2=", command2
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", shell_script_path]))
    scr.close()
    os.chmod(shell_script_path, 0777)

    return shell_script_path

# -----------------------------------------------------------------------

def submit_job(shell_script_path, my_type, options):
    if re.search("qsub", options.submission_command):
        qsub = "/sge/current/bin/lx-amd64/qsub"
        prog_name = "ep.merge." + str(my_type)
        command = qsub + " -V -N " + prog_name 
        command += " -A 093307 "
        if options.submission_command == "qsub":
            command += " -o /dev/null -e /dev/null " 
        command += shell_script_path
    else:
        command = "source " + shell_script_path
    print "shell_script_path=", shell_script_path
    print "command=", command
    os.system(command)

# -----------------------------------------------------------------------

if __name__ == "__main__":
   
    usage = "\nUsage: %prog input_file type [ shank ]"   

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = file_merging_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.compile:
        print "\nCompiling Matlab code ..."
        compile_code(options)

    if len(args) == 1:
        input_file_path = os.path.join(tonic_data, args[0])
        for file_type in ["ft", "ss"]:
            command = os.path.join(tonic_home, "HighPerformanceComputing", "TNC_HPC_MergeMatFiles.py")
            os.system(command + " " + input_file_path + " " + file_type)
    elif len(args) == 2 or len(args) == 3:   # wf or mc file
        input_file_path = os.path.join(tonic_data, args[0])
        prefix = input_file_path.split(".")[0]   
        file_type = args[1]
        shank = ""
        if len(args) == 3:
            shank = args[2]
        shell_script_path = create_shell_script(input_file_path, file_type, shank, options)
        submit_job(shell_script_path, file_type, options)
    elif options.compile:
        sys.exit(2)
    else:
        print usage
        sys.exit(2)
