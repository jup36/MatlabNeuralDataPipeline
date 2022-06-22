#! /usr/bin/python                        
#
# Copyright (C) 2013 by Howard Hughes Medical Institute.
#
# Purpose: extarct the solution with highest fmin 
#          from all files produced by cluster nodes

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

def extract_best_solution_command_line_parser(parser):
    parser.add_option("-c", "--compile",   action="store_true", dest="compile", help="compile Matlab code", metavar="compile", default=False)
    parser.add_option("-D", "--debug",     dest="debug", help="debugging mode; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-d", "--output_folder", dest="output_folder", help="folder where to place executable", metavar="output_folder", default=tonic_home + "/..")
    parser.add_option("-e", "--executable_name", dest="executable_name", help="name of matlab executable", metavar="executable_name", default="TNC_HPC_ExtractBestSolution")
    parser.add_option("-F", "--num_folds", dest="num_folds", help="# folds to be used in cross-validation", metavar="num_folds", default="2")
    parser.add_option("-M", "--num_components", dest="num_components", help="# of Gaussian components", metavar="num_components", default="")
    parser.add_option("-n", "--num_replicas", dest="num_replicas",  help="# replics of solution used to determine the best (def=3)", metavar="num_replicas", default="1")
    parser.add_option("-S", "--shunks", dest="shanks", help="shanks to be processed", metavar="num_shank", default="1,2,3,4,5,6,7,8")
    parser.add_option("-s", "--submission_command", dest="submission_command", help="source, qsub or qsub_debug", metavar="submission_command", default="qsub")
    parser.add_option("-T", "--segments",   dest="segments", help="time segment to be processed", metavar="segments", default="1,2,3,4,5,6,7,8")
    parser.add_option("-v", "--verbose",   action="store_true", dest="verbose",               help="increase the verbosity level of output", default=False)
    return parser

# -----------------------------------------------------------------------

def compile_code(options):
    input_path      = os.path.join(tonic_home,"HighPerformanceComputing", \
                                   "TNC_HPC_ExtractBestSolution.m")
    command = "mcc -m " + input_path + " -o " + options.executable_name + \
              " -d " + options.output_folder
    if options.verbose:
        command += " -v "
        print "command= ", command, "\n"
    os.system(command)
    os.system("chmod +x " + os.path.join(options.output_folder, \
                                         options.executable_name))

# -----------------------------------------------------------------------

def process_inputs(ft_file, node, dict_nodes_extract_best_solution, options):
    name_prefix     = ft_file[0:(len(ft_file)-7)]
    executable_path = os.path.join(tonic_home, "..", options.executable_name)

    node           = int(node)
    shank          = dict_nodes_extract_best_solution[node][0]
    segment        = dict_nodes_extract_best_solution[node][1]
    num_components = dict_nodes_extract_best_solution[node][2]
    fold           = dict_nodes_extract_best_solution[node][3]

    output_prefix = name_prefix + "_shank" + str(shank) \
                                + "_seg"   + str(segment) \
                                + "_mc"    + str(num_components) \
                                + "_gm_"   + str(fold)
    output_name = output_prefix + ".mat"
    
    command_ebs = executable_path
    command_rm  = "rm -f "
    for replica in range(1, int(options.num_replicas)+1):
        input_name   = output_prefix + "_" + str(replica) + ".mat"
        command_ebs += " " + input_name
        command_rm  += " " + input_name
    command_ebs     += " " + output_name                 

    if options.verbose:
        print "low level command_ebs=", command_ebs
    os.system(command_ebs)

    if not options.debug:
        print "In TNC_HPC_ExtractBestSolution.py: command_rm=", command_rm
        os.system(command_rm)

# ----------------------------------------------------------------------

def map_nodes_extract_best_solution(options):
    shanks   = parse_option(options.shanks)
    segments = parse_option(options.segments)
    M        = parse_option(options.num_components)
    nf       =          int(options.num_folds)

    dict_nodes_extract_best_solution = {}
    node = 0
    for sh in shanks:
        for sg in segments:
            for nc in M:
                for f in range(0, nf + 1):
                    # NOTE: f == 0 is interpreted as num_folds == 1
                    node += 1
                    dict_nodes_extract_best_solution[node] = [int(sh), int(sg), \
                                                              int(nc), int(f)]                                                         
    return (node, dict_nodes_extract_best_solution)

# ----------------------------------------------------------------------

def parse_option(option):
    if re.search("-", option) and re.search(",", option):
        sys.exit("Please, use either dash or comma in specifying " + option + "\n")
    if re.search("-", option):
        splt = option.split("-")
        items = range(int(splt[0]), int(splt[1])+1)
    else:
        if re.search(",", option):
            items  = option.split(",")
        else:
            items  = [option]
    return items

# -----------------------------------------------------------------------

if __name__ == "__main__":
   
    usage = "\nUsage: %prog ft_file node [options (-h to list)]\n"

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = extract_best_solution_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.compile:
        print "\nCompiling Matlab code ..."
        compile_code(options)
    if len(args) == 2:
        num_nodes, dict_nodes_extract_best_solution = \
            map_nodes_extract_best_solution(options)
        ft_file, node = args[0:2]
        if int(node) > num_nodes:
            print "Index of node=", node, " > num_nodes=", str(num_nodes), "  in TNC_HPC_ExtractBestSolution.py"
            sys.exit(2)
        process_inputs(ft_file, node, dict_nodes_extract_best_solution, options)
    elif options.compile:
        sys.exit(2)
    else:
        print usage
        sys.exit(2)
