#! /usr/bin/python                        
#
# Copyright (C) 2013 by Howard Hughes Medical Institute.
#

# ------------------------- imports -------------------------

import os
import shutil, commands
import sys, re, optparse

tonic_home = "/groups/zlatic/home/denisovg/work/Josh_Dudman/TONIC_v2" # default location for Sun Grid processing
tonic_data = "/groups/zlatic/home/denisovg/work/Josh_Dudman" 
if "TONIC_HOME" in os.environ.keys():
    tonic_home = os.environ['TONIC_HOME']
if "TONIC_DATA" in os.environ.keys():
    tonic_data = os.environ['TONIC_DATA']

# -----------------------------------------------------------------------

def cross_validation_command_line_parser(parser):
    parser.add_option("-c", "--compile",   action="store_true", dest="compile", help="compile Matlab code", metavar="compile", default=False)
    parser.add_option("-C", "--computational_method", dest="method", help="em (default), ms or kk", metavar="method", default="em")
    parser.add_option("-d", "--output_folder", dest="output_folder", help="folder where to place executable", metavar="output_folder", default=tonic_home + "/..")
    parser.add_option("-D", "--debug",     dest="debug", help="debugging mode; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-e", "--executable_name", dest="executable_name",  help="name of matlab executable", metavar="executable_name",  default="TNC_GM_CrossValidation")
    parser.add_option("-F", "--num_folds", dest="num_folds", help="# folds to be used in cross-validation", metavar="num_folds", default="2")
    parser.add_option("-i", "--ifold", dest="ifold", help="index of the fold to be tested",  metavar="ifold", default="")
    parser.add_option("-M", "--num_components", dest="num_components", help="#'s of Gaussian components considered in the models", metavar="num_components", default="")
    parser.add_option("-m", "--model", dest="model", help="model to be used (=1 or 2)",  metavar="model", default="1")
    parser.add_option("-S", "--shunks", dest="shanks", help="shanks to be processed", metavar="shanks", default="1,2,3,4,5,6,7,8")
    parser.add_option("-s", "--submission_command", dest="submission_command", help="source, qsub or qsub_debug", metavar="submission_command", default="qsub")
    parser.add_option("-T", "--time_segments",   dest="segments", help="time segment to be processed", metavar="segments", default="1,2,3,4,5,6,7,8")
    parser.add_option("-v", "--verbose",   action="store_true", dest="verbose",               help="increase the verbosity level of output", default=False)
    return parser

# -----------------------------------------------------------------------

def compile_code(options):
    input_path  = os.path.join(tonic_home,"GaussianMixtureModel", \
                               "TNC_GM_CrossValidation.m")
    command = "mcc -m "+ input_path + " -o " + options.executable_name + \
              " -d "    + options.output_folder
    if options.verbose:
        command += " -v "
        print "command= ", command, "\n"
    os.system(command)
    os.system("chmod +x " + os.path.join(options.output_folder, \
                                         options.executable_name))

# -----------------------------------------------------------------------

def process_inputs(ft_file_path, node, dict_nodes_cross_validation,\
                                     options):
    name_prefix = ft_file_path[0:(len(ft_file_path)-7)];
    executable_path = os.path.join(tonic_home, "..", options.executable_name)

    node           = int(node)
    shank          = dict_nodes_cross_validation[node][0]
    segment        = dict_nodes_cross_validation[node][1]
    num_components = dict_nodes_cross_validation[node][2]

    command_cv     = executable_path     + " " + ft_file_path + " " +\
                     str(shank)          + " " + str(segment) + " " +\
                     str(options.method) + " " + str(options.model) + " " +\
                     str(num_components) + " " + str(options.num_folds) 
    print "command_cv=", command_cv
    output_prefix = name_prefix + "_shank" + str(shank) \
                                + "_seg"   + str(segment) \
                                + "_mc"    + str(num_components) \
                                + "_gm"
    output_name = output_prefix + ".mat"
    command_rm  = "rm -f "
    for fold in range(1, int(options.num_folds)+1):
        input_name   = output_prefix + "_" + str(fold) + ".mat"
        command_rm  += " " + input_name
    command_cv      += " " + output_name
    if options.verbose:
        print "command_cv=", command_cv
    os.system(command_cv)
    if not options.debug:
        print "In TNC_HPC_CrossValidation.py: command_rm=", command_rm
        os.system(command_rm)

# ----------------------------------------------------------------------

def map_nodes_cross_validation(options):
    shanks   = parse_option(options.shanks)
    segments = parse_option(options.segments)
    M        = parse_option(options.num_components)

    dict_nodes_cross_validation = {}
    node = 0
    for sh in shanks:
        for sg in segments:
            for nc in M:
                node += 1
                dict_nodes_cross_validation[node] = [int(sh), int(sg), \
                                                     int(nc)]
    return (node, dict_nodes_cross_validation)

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
    parser = cross_validation_command_line_parser(parser)
    (options, args) = parser.parse_args()

    print "len(args)=", len(args), " agrs=", args
    print "options.submission_command=", options.submission_command
    if options.compile:
        print "\nCompiling Matlab code ..."
        compile_code(options)
    if len(args) == 2:           
        num_nodes, dict_nodes_cross_validation = \
            map_nodes_cross_validation(options)
        ft_file_path, node = args[0:2]
        print "node=", node, " num_nodes=", num_nodes
        if int(node) > num_nodes:
            print "Index of node > num_nodes in TNC_HPC_CrossValidation.py"
            sys.exit(2)
        process_inputs(ft_file_path, node, dict_nodes_cross_validation,\
                        options)
    elif options.compile:
        sys.exit(2)
    else:
        print usage
        sys.exit(2)
