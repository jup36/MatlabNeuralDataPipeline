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
    parser.add_option("-c", "--compile",   action="store_true", dest="compile", help="compile Matlab code", metavar="compile", default=False)
    parser.add_option("-C", "--computational_method", dest="method", help="em (default), ms or kk", metavar="method", default="em")
    parser.add_option("-d", "--output_folder", dest="output_folder", help="folder where to place executable", metavar="output_folder", default=tonic_home + "/..")
    parser.add_option("-e", "--executable_name", dest="executable_name", help="executable", metavar="executable_name", default="TNC_GM_SummarizeResults")
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
    input_path      = os.path.join(tonic_home,"GaussianMixtureModel", \
                                   "TNC_GM_SummarizeResults.m")
    command = "mcc -m " +\
              input_path + " -o " + options.executable_name + \
              " -d " + options.output_folder
    if options.verbose:
        command += " -v "
        print "command= ", command, "\n"
    os.system(command)
    os.system("chmod +x " + os.path.join(options.output_folder, \
                                         options.executable_name))

# -----------------------------------------------------------------------------

def summarize_results(input_files, options):
    command = os.path.join(tonic_data, options.executable_name)
    for input_file in input_files:
        command += " " + input_file
    if options.verbose:
        print "command=", command
    os.system(command)

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
#       print "\nargs=", args, "\n"
        input_files = []
        for arg in args:
            input_files1 = parse_input_data(arg, tonic_data)
            for if1 in input_files1:
                if not if1 in input_files:
                    input_files.append(if1)
        if options.verbose:
            print "input files=", input_files
        summarize_results(input_files, options)
    elif options.compile:
        sys.exit(2)
    else:
        parser.print_usage()
        sys.exit(2)
