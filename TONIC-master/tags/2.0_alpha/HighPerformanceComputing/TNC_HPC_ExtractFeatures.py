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

def feature_extraction_command_line_parser(parser):
    parser.add_option("-A", "--project",   dest="project_code", help="project code to be used with qsub", metavar="project_code", default="093307")
    parser.add_option("-c", "--compile",   action="store_true", dest="compile", help="compile Matlab code", metavar="compile", default=False)
    parser.add_option("-e", "--executable_name", dest="executable_name", help="executable", metavar="executable_name", default="TNC_SSP_ExtractFeatures")
    parser.add_option("-d", "--output_folder", dest="output_folder", help="folder containing executable", metavar="output_folder", default=tonic_home + "/..")
    parser.add_option("-D", "--debug",     dest="debug", help="debugging mode; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-f", "--fake_data", action="store_true", dest="fake_data", help="specify that fake data is processed", default=False)
    parser.add_option("-M", "--scm",   action="store_true", dest="scm",   help="subtract common mean when processing traces", default=False)
    parser.add_option("-n", "--num_electrodes", dest="num_electrodes", help="# electrodes used to generate data", metavar="num_electrodes", default="")
    parser.add_option("-r", "--snr",       dest="snr", help="threshold value of signal-to-noise ratio", metavar="snr", default="");
    parser.add_option("-s", "--sub_command", dest="submission_command", help="source, qsub or qsub_debug", metavar="submission_command", default="qsub")
    parser.add_option("-t", "--array_type",dest="array_type",  help="electrode array type", metavar="which shank to process", default="")
    parser.add_option("-T", "--num_segments",   dest="num_segments", help="# of time segment to be ptocessed", metavar="num_segments", default="8") 
    parser.add_option("-v", "--verbose",   action="store_true", dest="verbose",   help="increase the verbosity level of output", default=False)
    return parser


# -----------------------------------------------------------------------

# Takes comma-separated list of files, possibly containing wild cards
# All the files are assumed to be located in folder tonic_data
# Outputs a list of real names of the files

def parse_input_data(input_item, tonic_data):
    wild_cards = input_item.split(",")

    # Identify data dir
    input_files = []

    for wc in wild_cards:
        if not re.search("\*", wc): 
            my_path = os.path.join(tonic_data, wc)
            if os.path.exists(my_path):
                input_files.append(wc)
        else:    
            fragments = wc.split("*")
            for file in os.listdir(tonic_data):
                file_matched = True
                f = fragments[len(fragments)-1]
#               print "fragments=", fragments, " file=", file, " beg=", file[0:len(fragments[0])], " end=", file[(len(file)-len(f)):len(file)]
                for f in fragments:
                    if not re.search(f, file):
                        file_matched = False
                        break
                    if fragments.index(f) == 0:
                        if not file[0:len(fragments[0])] == f:
                           file_matched = False
                           break
                    if fragments.index(f) == len(fragments)-1:
                        if not file[(len(file)-len(fragments[len(fragments)-1])):len(file)] == f:
                           file_matched = False
                           break
                if file_matched:
                    input_files.append(file)
                
    return input_files

# -----------------------------------------------------------------------

def compile_code(options):
    input_path      = os.path.join(tonic_home,"EventDetectionClassification", \
                                   "TNC_SSP_ExtractFeatures.m")
    command = "mcc -m -N -p signal -p optim -p stats -R -singleCompThread " +\
              input_path + " -o " + options.executable_name + \
              " -d " + options.output_folder
    if options.verbose:
        command += " -v "
        print "command= ", command, "\n"
    os.system(command)
    os.system("chmod +x " + os.path.join(options.output_folder, \
                                         options.executable_name))

# -----------------------------------------------------------------------

def create_high_level_shell_script(outfolderpath, input_file_path, options):
    shell_script_path = os.path.join(outfolderpath, "Shell_script_level2.sh")
    scr = open(shell_script_path, 'wt')
    scr.write("#!/usr/bash\n")
    for i in range(1, int(options.num_segments)+1):
        shell_script_command = os.path.join(tonic_home, "HighPerformanceComputing", \
                                                        "TNC_HPC_ExtractFeatures.py")
        shell_script_command += " " + input_file_path[0]\
                             +  " " + str(i)\
                             +  " -A " + options.project_code \
                             +  " -s " + options.submission_command \
                             +  " -T " + options.num_segments
        if options.verbose:
            shell_script_command += " -v "
        if options.scm:
            shell_script_command += " -M "
        shell_script_command += "\n"
        scr.write(shell_script_command)
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", shell_script_path])) 
    scr.close()
    os.chmod(shell_script_path, 0777)
    return shell_script_path

# ----------------------------------------------------------------------

def create_array_job_script(outfolderpath, input_file_path, options):
    tokens = input_file_path.split("/")
    input_file = tokens[len(tokens)-1]
    shell_script_path = os.path.join(outfolderpath, \
                        "Shell_script_array_job." + input_file + ".sh")
    executable_path = os.path.join(tonic_home, "..", options.executable_name)
    if options.array_type == "NN_b64":
        num_shanks = 8
    elif options.array_type == "tetrode":
        num_shanks = 1
    else:
        num_shanks = int(options.num_electrodes)
    num_nodes = num_shanks * int(options.num_segments) 
    print "num_nodes=", num_nodes

    scr = open(shell_script_path, 'wt')
    scr.write("#!/usr/bash\n")
#   scr.write("export MCR_CACHE_ROOT=/tmp/$USER/mcr_cache_root.$JOB_ID.$SGE_TASK_ID\n")
#   scr.write("mkdir -p $MCR_CACHE_ROOT\n")
    if re.search("qsub", options.submission_command):
        scr.write("#$ -t  1-" + str(num_nodes) + "\n")
        shell_scr_command =  executable_path + " '" + input_file_path + "' '" \
            + options.array_type + "' " + options.num_segments + " $SGE_TASK_ID " \
            + " fakeData " + str(int(options.fake_data)) \
            + " snr "      + str(options.snr) \
            + " scm "      + str(options.scm) \
            + " verbose " + str(int(options.verbose)) + " \n"
        print "shell_scr_command=", shell_scr_command
        scr.write(shell_scr_command)
    else:
        for i in range(0, num_nodes):
            shell_scr_command =  executable_path + " '" + input_file_path + "' '" \
            + options.array_type + "' " + options.num_segments + " " + str(i+1) \
            + " fakeData " + str(int(options.fake_data)) \
            + " snr "      + str(options.snr) \
            + " scm "      + str(options.scm) \
            + " verbose " + str(int(options.verbose)) + " \n"
            scr.write(shell_scr_command)
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", shell_script_path])) 
    scr.close()
    os.chmod(shell_script_path, 0777)
    return shell_script_path

# -----------------------------------------------------------------------------

def create_epilog_script(outfolderpath, input_file_path, options):
    tokens = input_file_path.split("/")
    input_file = tokens[len(tokens)-1]
    epilog_script_path = os.path.join(outfolderpath, \
                                      "Epilog_script_feature_extraction" +\
                                      input_file + ".sh")

    scr = open(epilog_script_path, 'wt')
    scr.write("#!/usr/bash\n")
    prefix = input_file.split(".")[0]
    for my_type in ["ft", "ss"]:
        command = os.path.join(tonic_home, "HighPerformanceComputing", \
                               "TNC_HPC_MergeMatFiles.py") + " "
        command += input_file_path + " " + my_type 
        command += " -s source "
        if options.verbose:
            command += " -v "
        if options.debug:
            command += " -D "
        scr.write(command + "\n")
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", epilog_script_path]))
    scr.close()
    return epilog_script_path

#------------------------------------------------------------------------------

def submit_array_job(shell_script_path, epilog_script_path, options):
    qsub = "/sge/current/bin/lx-amd64/qsub"
    prog_name = "xf.array"
    if re.search("qsub", options.submission_command):
        command = qsub + " -V -N " + prog_name
        command += " -pe batch 9 "
        if len(options.project_code) > 0:
            command += "  -A " + options.project_code
        if options.submission_command == "qsub":
            command += " -o /dev/null -e /dev/null "
        command += " " + shell_script_path
        res = commands.getstatusoutput(command)
        array_jobid = (res[1].split()[2]).split(".")[0]
        ep_command = qsub + " -V -N xf.epilog -hold_jid " + array_jobid 
        if len(options.project_code) > 0:
            ep_command += " -A " + options.project_code
        if options.submission_command == "qsub":
            ep_command += " -o /dev/null -e /dev/null " 
        ep_command += " " + epilog_script_path
        os.system(ep_command)
    else:
        command = "source  "
        print "command =", command 
        os.system(command)
        ep_command = "source " + epilog_script_path
        os.system(ep_command)
    if options.verbose:
        print "array job command=", command, "\nres=", res
        print "ep_command=", ep_command

    # Submit epilog script
    return array_jobid

# -----------------------------------------------------------------------------

def extract_features_high_level(input_file_path, options):
    outfolderpath = tonic_data
    shell_script_path = \
        create_array_job_script(outfolderpath, input_file_path, options)
    epilog_script_path = \
        create_epilog_script(outfolderpath, \
                             input_file_path, options)
    if options.verbose:
        print "array_job_script_path=", shell_script_path
        print "epilog_script_path=", epilog_script_path
    submit_array_job(shell_script_path, epilog_script_path, options)

# ----------------------------------------------------------------------

if __name__ == "__main__":
   
    usage = "Usage: \n\
    %prog input_file(s) [options (-h to list)]"   

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = feature_extraction_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.compile:
        print "\nCompiling Matlab code ..."
        compile_code(options)
    if len(args) > 0:
        # Reset options
        if len(options.snr) == 0:
            options.snr = "-1";
        if len(options.array_type) == 0:
            sys.exit("Please, specify an array type with option -a")
        elif options.array_type == "NN_b64":
            options.num_electrodes = "64"
        elif options.array_type == "tetrode":
            options.num_electrodes = "4"
        else:
            if len(options.num_electrodes) == 0:
                sys.exit("Please, specify a number of electrodes with option -n")
        input_files = []
        for arg in args:
            input_files1 = parse_input_data(arg, tonic_data)
            for if1 in input_files1:
                if not if1 in input_files:
                    input_files.append(if1)
        print "input files=", input_files
        for input_file in input_files:
            input_file_path = os.path.join(tonic_data, input_file)
            extract_features_high_level(input_file_path, options)
    elif options.compile:
        sys.exit(2)
    else:
        parser.print_usage()
        sys.exit(2)
