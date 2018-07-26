#! /usr/local/python-2.7.6/bin/python     
#
# Copyright (C) 2013 by Howard Hughes Medical Institute.
#

# Purpose: submit Matlab spike sorting job(s) to cluster

# ------------------------- imports -------------------------

import os, glob, gc
import shutil, commands
import sys, re, optparse
import scipy.io

tonic_home = os.environ['TONIC_HOME']
tonic_data = os.environ['TONIC_DATA']

# ----------------------------------------------------------------------

def spike_sorting_command_line_parser(parser):
    parser.add_option("-A","--project",   dest="project_code", help="project code to be used with qsub", metavar="project_code", default="093307")
    parser.add_option("-a","--array_type",dest="array_type",  help="electrode array type", metavar="which shank to process", default="NN_b64")
    parser.add_option("-c","--compile",   action="store_true", dest="compile", help="compile Matlab code", metavar="compile", default=False)
    parser.add_option("-C","--compile_all", action="store_true", dest="compile_all", help="compile all spike sorting Matlab code", metavar="compilei_all", default=False)
    parser.add_option("-t","--computational_tool", dest="method", help="em (default), ms or kk", metavar="method", default="em")
    parser.add_option("-D","--debug",     dest="debug", help="debugging mode; don't delete shell scripts", action="store_true", default=False)
    parser.add_option("-d","--output_folder", dest="output_folder", help="folder where to place executable", metavar="output_folder", default=tonic_home + "/..")
    parser.add_option("-e","--executable_name", dest="executable_name", help="name of matlab executable", metavar="executable_name", default="TNC_GM_SortSpikes")
    parser.add_option("-f","--fake_data", action="store_true", dest="fake_data", help="specify that fake data is processed", default=False)
    parser.add_option("-F","--num_folds", dest="num_folds", help="# folds to be used in cross-validation", metavar="num_folds", default="3")
    parser.add_option("-g","--guess_method", dest="guess_method",  help="method to compute initial guess for solution", metavar="init_guess", default="km")
    parser.add_option("-i","--ifold", dest="ifold", help="index of the fold to be tested",  metavar="ifold", default="")
    parser.add_option("-I","--max_iter", dest="max_iter", help="max. number of iterations by spike sorting algorithm",metavar="max_iter", default="100")
    parser.add_option("-m","--model", dest="model", help="model to be used (=1 or 2)",  metavar="model", default="1")
    parser.add_option("-M","--num_components", dest="num_components", help="comma- or dash- separated range of considered Gaussian components, default='1,2,3,4,5,6,7'", metavar="num_components", default="1,2,3,4,5,6,7")
    parser.add_option("-n","--num_replicas", dest="num_replicas",  help="# replicas of solution used to determine the best, default=1", metavar="num_replicas", default="1")
    parser.add_option("-N","--inode",   dest="inode", help="id of the cluster node to be used", metavar="inode",  default="")
    parser.add_option("-p", "--processing_start",dest="processing_start",help="start from rem(1),ebs(2),cvl(3),bnc(4),or epi(5), default=5 ",metavar="processing_start",default=1)
    parser.add_option("-P", "--processing_end",dest="processing_end",help="end with rem(1),ebs(2),cvl(3),bnc(4),or epi(5), default=5 ",metavar="processing_end",default=5)
    parser.add_option("-S","--shunks", dest="shanks", help="comma- or dash-separated list of shanks to be processed, default=all", metavar="shanks", default="all")
    parser.add_option("-s","--submission_command", dest="submission_command", help="source, qsub or qsub_debug, default=qsub", metavar="submission_command", default="qsub")
    parser.add_option("-T","--segments",dest="segments",help="comma- or dash-separated list of time segment(s) to be processed, default=all",metavar="segments",default="all")
    parser.add_option("-v","--verbose",   action="store_true", dest="verbose", help="increase the verbosity level of output", default=False)
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
    print "Compiling TNC_GM_SortSpikes.m ..."
    input_path = os.path.join(tonic_home,"GaussianMixtureModel", 
                              "TNC_GM_SortSpikes.m")
    command = "mcc -m -N -p globaloptim -p optim -p stats -R -singleCompThread "\
              + input_path + " -o " + options.executable_name + \
              " -d " + os.path.join(tonic_home, "Bin")
    if options.verbose:
        command += " -v "
        print "command= ", command, "\n"
    os.system(command)
    os.system("chmod +x " + os.path.join(tonic_home, "Bin", options.executable_name))
    if options.compile_all:
        print "Compiling TNC_HPC_ExtractBestSolution.m ..."
        input_path = os.path.join(tonic_home,"HighPerformanceComputing", \
                                  "TNC_HPC_ExtractBestSolution.m")
        command = "mcc -m " + input_path + " -o " + "TNC_HPC_ExtractBestSolution"  + \
                  " -d " + os.path.join(tonic_home, "Bin")
        if options.verbose:
            command += " -v "
            print "command= ", command, "\n"
        os.system(command)
        os.system("chmod +x " + os.path.join(tonic_home, "Bin", "TNC_HPC_ExtractBestSolution"))
        # --------------------------------------------------------------
        print "Compiling TNC_GM_BestNumComponents.m ..."
        input_path = os.path.join(tonic_home,"GaussianMixtureModel", \
                                  "TNC_GM_BestNumComponents.m")
        command = "mcc -m " + input_path + " -o " + "TNC_GM_BestNumComponents"  + \
                  " -d " + os.path.join(tonic_home, "Bin")
        if options.verbose:
            command += " -v "
            print "command= ", command, "\n"
        os.system(command)
        os.system("chmod +x " + os.path.join(tonic_home, "Bin", "TNC_GM_BestNumComponents"))
        # --------------------------------------------------------------
        print "Compiling TNC_GM_CrossValidation.m ..."
        input_path = os.path.join(tonic_home,"GaussianMixtureModel", \
                                  "TNC_GM_CrossValidation.m")
        command = "mcc -m " + input_path + " -o " + "TNC_GM_CrossValidation"  + \
                  " -d " + os.path.join(tonic_home, "Bin")
        if options.verbose:
            command += " -v "
            print "command= ", command, "\n"
        os.system(command)
        os.system("chmod +x " + os.path.join(tonic_home, "Bin", "TNC_GM_CrossValidation"))
        # --------------------------------------------------------------
        print "\nCompiling TNC_HPC_MergeMatFiles.m ..."
        input_path  = os.path.join(tonic_home,"HighPerformanceComputing", \
                                   "TNC_HPC_MergeMatFiles.m")
        command = "mcc -m " + input_path + " -o " + "TNC_HPC_MergeMatFiles" + \
                  " -d " + os.path.join(tonic_home, "Bin")
        if options.verbose:
            command += " -v "
            print "command= ", command, "\n"
        os.system(command)
        os.system("chmod +x " + os.path.join(tonic_home, "Bin", "TNC_HPC_MergeMatFiles")) 
        
# -----------------------------------------------------------------------

def create_rem_script(ft_file_path, outfolderpath, options):
    tokens = ft_file_path.split("/")
    input_file = tokens[len(tokens)-1]
    run_em_script_path = os.path.join(outfolderpath, \
                        "Run_em." + input_file + ".sh")
    executable_path = os.path.join(tonic_home, "HighPerformanceComputing", \
                                   "TNC_HPC_SortSpikes.py")
    scr = open(run_em_script_path, 'wt')
    num_nodes_rem, dict_nodes_rem = map_nodes_rem(ft_file_path, options)
    if options.verbose:
        print "num_nodes_rem=", num_nodes_rem

    scr.write("#!/usr/bash\n")
    scr.write("#$ -t 1-" + str(num_nodes_rem) + "\n")
    my_segments_str = str(options.segments)
    my_shanks_str   = str(options.shanks)
    if options.segments == "all":
        num_segments = get_num_segments_from_data(ft_file_path)
        my_segments_str = ",".join([str(i) for i in range(1,num_segments+1)])
    if options.shanks == "all":
        num_shanks = get_num_shanks_from_data(ft_file_path)
        my_shanks_str = ",".join([str(i) for i in range(1,num_shanks+1)])
    command_rem = executable_path + " '" + ft_file_path + "' " \
                   + " -N $SGE_TASK_ID "\
                   + " -t " + options.method \
                   + " -T " + my_segments_str  \
                   + " -S " + my_shanks_str \
                   + " -M " + str(options.num_components)\
                   + " -n " + str(options.num_replicas)\
                   + " -m " + str(options.model) \
                   + " -I " + str(options.max_iter)\
                   + " -F " + str(options.num_folds)
    if options.verbose:
        command_rem += " -v"
    if options.fake_data:
        command_rem += " -f "
    command_rem += " -e " + options.executable_name
    command_rem += "\n"
    if options.verbose:
        print "command_rem=", command_rem
    scr.write(command_rem)
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", run_em_script_path])) 
    scr.close()
    os.chmod(run_em_script_path, 0777)
    return run_em_script_path

# -----------------------------------------------------------------------------

def create_extract_best_solution_script(ft_file_path, outfolderpath, options):
    tokens = ft_file_path.split("/")
    input_file = tokens[len(tokens)-1]
    extract_best_solution_script_path = os.path.join(outfolderpath, \
                         "Extract_best_solution_script." + input_file + ".sh")
    executable_path = os.path.join(tonic_home, "HighPerformanceComputing", \
                                   "TNC_HPC_ExtractBestSolution.py")
    scr = open(extract_best_solution_script_path, 'wt')
    num_nodes_extract_best_solution, dict_nodes_extract_best_solution = \
        map_nodes_extract_best_solution(ft_file_path, options)
    scr.write("#!/usr/bash\n")
    scr.write("#$ -t 1-" + str(num_nodes_extract_best_solution) + "\n")
    my_segments_str = str(options.segments)
    my_shanks_str   = str(options.shanks)
    if options.segments == "all":
        num_segments = get_num_segments_from_data(ft_file_path)
        my_segments_str = ",".join([str(i) for i in range(1,num_segments+1)])
    if options.shanks == "all":
        num_shanks = get_num_shanks_from_data(ft_file_path)
        my_shanks_str = ",".join([str(i) for i in range(1,num_shanks+1)])
    command_ebs = executable_path + " '" + ft_file_path + "' " \
                           " $SGE_TASK_ID " + \
                           " -T " + my_segments_str + \
                           " -S " + my_shanks_str   + \
                           " -M " + str(options.num_components) +\
                           " -n " + options.num_replicas + \
                           " -s " + options.submission_command + \
                           " -F " + options.num_folds
    if options.debug:
        command_ebs += " -D "
    if options.verbose:
        command_ebs += " -v "
    if options.verbose:
        print "In create_extract_best_solution_script: command_ebs=", command_ebs
    scr.write(command_ebs + "\n")
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", extract_best_solution_script_path]))
    scr.close()
    return extract_best_solution_script_path

# -----------------------------------------------------------------------------

def create_cross_validation_script(ft_file_path, outfolderpath, options):
    tokens = ft_file_path.split("/")
    input_file = tokens[len(tokens)-1]
    cross_validation_script_path = os.path.join(outfolderpath, \
                         "Cross_validation_script." + input_file + ".sh")
    executable_path = os.path.join(tonic_home, "HighPerformanceComputing", \
                                   "TNC_HPC_CrossValidation.py")              
    scr = open(cross_validation_script_path, 'wt')
    num_nodes_cross_validation, dict_nodes_cross_validation = \
        map_nodes_cross_validation(ft_file_path, options)
    scr.write("#!/usr/bash\n")
    scr.write("#$ -t 1-" + str(num_nodes_cross_validation) + "\n")

    my_segments_str = str(options.segments)
    my_shanks_str   = str(options.shanks)
    if options.segments == "all":
        num_segments = get_num_segments_from_data(ft_file_path)
        my_segments_str = ",".join([str(i) for i in range(1,num_segments+1)])
    if options.shanks == "all":
        num_shanks = get_num_shanks_from_data(ft_file_path)
        my_shanks_str = ",".join([str(i) for i in range(1,num_shanks+1)])

    command_cv  = executable_path + " " + ft_file_path +\
                            " $SGE_TASK_ID " + \
                            " -F " + options.num_folds + \
                            " -T " + my_segments_str +\
                            " -S " + my_shanks_str   +\
                            " -s " + options.submission_command + \
                            " -M " + str(options.num_components) +\
                            " -m " + options.model 
    if options.debug:
        command_cv  += " -D "
    if options.verbose:
        command_cv  += " -v "
        print "command_cv=", command_cv
    scr.write(command_cv  + "\n")
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", cross_validation_script_path]))
    scr.close()
    return cross_validation_script_path

#------------------------------------------------------------------------------

def create_best_num_components_script(ft_file_path, outfolderpath, options):
    tokens = ft_file_path.split("/")
    input_file = tokens[len(tokens)-1]
    best_num_components_script_path = os.path.join(outfolderpath, \
                         "Best_num_components_script." + input_file + ".sh")
    executable_path = os.path.join(tonic_home, "HighPerformanceComputing", \
                                   "TNC_HPC_BestNumComponents.py")
    scr = open(best_num_components_script_path, 'wt')
    num_nodes_best_num_components, dict_nodes_best_num_components = \
        map_nodes_best_num_components(ft_file_path, options)
    scr.write("#!/usr/bash\n")
    scr.write("#$ -t 1-" + str(num_nodes_best_num_components) + "\n")

    my_segments_str = str(options.segments)
    my_shanks_str   = str(options.shanks)
    if options.segments == "all":
        num_segments = get_num_segments_from_data(ft_file_path)
        my_segments_str = ",".join([str(i) for i in range(1,num_segments+1)])
    if options.shanks == "all":
        num_shanks = get_num_shanks_from_data(ft_file_path)
        my_shanks_str = ",".join([str(i) for i in range(1,num_shanks+1)])

    command_bnc = executable_path + " " + ft_file_path +\
                            " $SGE_TASK_ID " + \
                            " -t " + options.method +\
                            " -F " + options.num_folds + \
                            " -M " + options.num_components +\
                            " -m " + str(options.model) +\
                            " -T " + my_segments_str +\
                            " -S " + my_shanks_str   +\
                            " -s " + options.submission_command + \
                            " -F " + str(options.num_folds)
    if options.debug:
        command_bnc += " -D "
    if options.verbose:
        command_bnc += " -v "
        print "command_bnc=", command_bnc
    scr.write(command_bnc + "\n")
    if not options.debug:
        scr.write("%s '%s' \n" % tuple(["rm -f ", best_num_components_script_path]))
    scr.close()
    return best_num_components_script_path

#------------------------------------------------------------------------------

def create_epilog_script(ft_file_path, outfolderpath, options):
    tokens = ft_file_path.split("/")
    ft_file = tokens[len(tokens)-1]
    name_prefix = ft_file[0:(len(ft_file)-7)];
    if not options.debug:
        output_file = os.path.join(outfolderpath, name_prefix + "_gm.mat")
    else:
        output_file = os.path.join(outfolderpath, name_prefix + \
                               "_nf" + str(options.num_folds) + "_gm.mat")
    epilog_script_path = os.path.join(outfolderpath, \
                         "Epilog_merge_mat_files." + ft_file + ".sh")
    scr = open(epilog_script_path, 'wt')
    scr.write("#!/usr/bash\n")
    command_mmf = os.path.join(tonic_home, "Bin", "TNC_HPC_MergeMatFiles")
    command_rm  = "rm -f "
    command_rm2 = "rm -f "

    if options.shanks == "all":
        num_shanks = get_num_shanks_from_data(ft_file_path)
        shanks = range(1, num_shanks+1)
    else:
        shanks   = parse_option(options.shanks)
    if options.segments == "all":
        num_segments = get_num_segments_from_data(ft_file_path)
        segments = range(1, num_segments+1)
    else:
        segments = parse_option(options.segments)

    for shank in shanks:
        for seg in segments:
            input_file = os.path.join(outfolderpath, name_prefix + \
                                      "_shank" + str(shank) + \
                                      "_seg"   + str(seg) + "_gm.mat")
            command_mmf += " " + input_file
            command_rm  += " " + input_file
    command_mmf += " " + output_file + " " + str(int(options.debug))
    scr.write(command_mmf + "\n")
    if not options.debug:
        scr.write(command_rm  + "\n")
        scr.write("%s '%s' \n" % tuple(["rm -f ", epilog_script_path]))
    scr.close()
    return epilog_script_path

#------------------------------------------------------------------------------

def submit_rem_job(run_em_script_path, command, options):
    prog_name1 = "ss.rem"

    if options.submission_command == "qsub":
        command1 = command + " -V -N " + prog_name1 + " -o /dev/null -e /dev/null "
    else:
        command1 = command + " -V -N " + prog_name1
    # Note: option -hold_jid <job_id> should precede the script name!!!
    #      (otherwise, this option will be ignored)
    command1 += " " + run_em_script_path
    res1 = commands.getstatusoutput(command1)
    jobid1 = (res1[1].split()[2]).split(".")[0]

    print "\nSubmit rum em command=", command1, "\nres1=", res1
    print "jobid1=", jobid1
    return jobid1

#------------------------------------------------------------------------------

def submit_ebs_job(ebs_script_path, command, jobid1, options):
    prog_name2 = "ss.ebs"
    if options.submission_command == "qsub":
        command2 = command + " -V -N " + prog_name2 +\
                             " -o /dev/null -e /dev/null "
    else:
        command2 = command + " -V -N " + prog_name2
    if jobid1 > 0:
        command2 += " -hold_jid " + jobid1
    command2 += " " +  ebs_script_path
    res2 = commands.getstatusoutput(command2)
    jobid2 = (res2[1].split()[2]).split(".")[0]

    print "\nSubmit ebs job command=", command2, "\nres2=", res2
    print "jobid2=", jobid2
    return jobid2

#------------------------------------------------------------------------------

def submit_cvl_job(cvl_script_path, command, jobid2, options):
    prog_name3 = "ss.cvl"
    if options.submission_command == "qsub":
        command3 = command + " -V -N " + prog_name3 +\
                             " -o /dev/null -e /dev/null "
    else:
        command3 = command + " -V -N " + prog_name3
    if jobid2 > 0:
        command3 += " -hold_jid " + jobid2
    command3 += " " +  cvl_script_path
    res3 = commands.getstatusoutput(command3)
    jobid3 = (res3[1].split()[2]).split(".")[0]

    print "\nSubmit cvl job command=", command3, "\nres3=", res3
    print "jobid3=", jobid3
    return jobid3

#------------------------------------------------------------------------------

def submit_bnc_job(bnc_script_path, command, jobid3, options):
    prog_name4 = "ss.bnc"
    if options.submission_command == "qsub":
        command4 = command + " -V -N " + prog_name4 +\
                             " -o /dev/null -e /dev/null "
    else:
        command4 = command + " -V -N " + prog_name4
    if jobid3 > 0:
        command4 += " -hold_jid " + jobid3
    command4 += " " +  bnc_script_path
    res4 = commands.getstatusoutput(command4)
    jobid4 = (res4[1].split()[2]).split(".")[0]

    print "\nSubmit bnc job command=", command4, "\nres4=", res4
    print "jobid4=", jobid4
    return jobid4

#------------------------------------------------------------------------------

def submit_epilog_job(epilog_script_path, command, jobid4, options):
    prog_name5 = "ss.epilog"
    if options.submission_command == "qsub":
        command5 = command + " -V -N " + prog_name5 +\
                             " -o /dev/null -e /dev/null "
    else:
        command5 = command + " -V -N " + prog_name5
    if jobid4 > 0:
        command5 += " -hold_jid " + jobid4
    command5 += " " +  epilog_script_path
    os.system(command5)

    if options.verbose:
        print "\nSubmit epilog job command=", command5

#------------------------------------------------------------------------------

def submit_array_jobs(rem_script_path, \
                      ebs_script_path,\
                      cvl_script_path,\
                      bnc_script_path,\
                      epilog_script_path, options):
    qsub = "/sge/current/bin/lx-amd64/qsub"
    if re.search("qsub", options.submission_command):
        base_command = qsub
    else:
        base_command = options.submission_command
    if len(options.project_code) > 0:
        base_command += "  -A " + options.project_code

    prog_name5 = "ss.epilog"

    jobid1 = 0
    if int(options.processing_start) == 1 and int(options.processing_end) >= 1:
        jobid1 = submit_rem_job(rem_script_path, base_command, options)
        if options.method == "kk":
            return jobid1

    jobid2 = 0
    if int(options.processing_start) <= 2 and int(options.processing_end) >= 2:
        jobid2 = submit_ebs_job(ebs_script_path, base_command, jobid1, options)

    jobid3 = 0
    if int(options.processing_start) <= 3 and int(options.processing_end) >= 3:
        jobid3 = submit_cvl_job(cvl_script_path, base_command, jobid2, options)

    jobid4 = 0
    if int(options.processing_start) <= 4 and int(options.processing_end) >= 4:
        jobid4 = submit_bnc_job(bnc_script_path, base_command, jobid3, options)

    if int(options.processing_start) <= 5:
        submit_epilog_job(epilog_script_path,    base_command, jobid4, options)

    return 

# -----------------------------------------------------------------------------

def sort_spikes_high_level(ft_file_path, method, num_nodes_rem, options):
    outfolderpath = tonic_data
    run_em_script_path = ""
    extract_best_solution_script_path = ""
    cross_validation_script_path      = ""
    best_num_components_script_path   = ""
    epilog_script_path                = ""
    if int(options.processing_start) <= 1 and int(options.processing_end) >= 1:
        run_em_script_path = \
            create_rem_script(ft_file_path, outfolderpath, options)
        if options.verbose:
            print "run_em_script_path=", run_em_script_path
    if not method == "kk":
        if int(options.processing_start) <= 2 and int(options.processing_end) >= 2:
            extract_best_solution_script_path = \
                create_extract_best_solution_script(ft_file_path, outfolderpath, options)
            if options.verbose:
                print "extract_best_solution_script_path=", extract_best_solution_script_path
        if int(options.processing_start) <= 3 and int(options.processing_end) >= 3:
            cross_validation_script_path = \
                create_cross_validation_script(ft_file_path, outfolderpath, options)
            if options.verbose:
                print "cross_validation_script_path=", cross_validation_script_path
        if int(options.processing_start) <= 4 and int(options.processing_end) >= 4:
            best_num_components_script_path = \
                create_best_num_components_script(ft_file_path, outfolderpath, options)
            if options.verbose:
                print "best_num_components_script_path=", best_num_components_script_path
        if int(options.processing_start) <= 5 and int(options.processing_end) >= 5:
            epilog_script_path = \
                create_epilog_script(ft_file_path, outfolderpath, options)
            if options.verbose:
                print "epilog_script_path=", epilog_script_path
    submit_array_jobs(run_em_script_path, \
                      extract_best_solution_script_path,\
                      cross_validation_script_path, \
                      best_num_components_script_path,\
                      epilog_script_path, options)

# -----------------------------------------------------------------------------

# Run spike sorting code for a given node

def rem_low_level(ft_file_path, node, dict_nodes_rem, options):
    outfolderpath = tonic_data
    executable_path = os.path.join(tonic_home, "Bin", options.executable_name)
    M = [ int(m) for m in parse_option(options.num_components) ]
    shank          = dict_nodes_rem[node][0]
    segment        = dict_nodes_rem[node][1]
    num_components = dict_nodes_rem[node][2]
    fold           = dict_nodes_rem[node][3]
    replica        = dict_nodes_rem[node][4]
    command = executable_path     + " " + ft_file_path + " " + \
              str(segment)        + " " + str(shank)   + " " + \
              str(num_components) + " " + str(replica) + " " + \
              " model "        + str(options.model)    + \
              " guess_method " + str(options.guess_method)  + \
              " verbose "      + str(int(options.verbose))  + \
              " maxIter "      + str(int(options.max_iter)) + \
              " iFold "        + str(fold) + \
              " numFolds "     + str(int(options.num_folds))
    if options.fake_data:
        command += " fakeData  1 " 
    if options.verbose:
        print "low level command_rem=", command
    os.system(command)

# ----------------------------------------------------------------------

def parse_option(option):
    items = []
    if re.search(",", option):
        split1 = option.split(",")
        for item1 in split1:
            if re.search("-", item1):
                bounds = item1.split("-")
                split2 = range(int(bounds[0]),int(bounds[1])+1)
                for item2 in split2:
                    if int(item2) not in items:
                        items.append(int(item2))
            else:
                if int(item1) not in items:
                    items.append(int(item1))
    elif re.search("-", option):
        bounds = option.split("-")
        split1 = range(int(bounds[0]),int(bounds[1])+1)
        for item1 in split1:
            if int(item1) not in items:
                items.append(int(item1))
    else:
        items  = [int(option)]
    return items   

# ----------------------------------------------------------------------

def map_nodes_rem(ft_file_path, options):
    if options.shanks == "all":
        num_shanks = get_num_shanks_from_data(ft_file_path)
        shanks = range(1, num_shanks+1)
    else:
        shanks   = parse_option(options.shanks)
    if options.segments == "all":
        num_segments = get_num_segments_from_data(ft_file_path)
        segments = range(1, num_segments+1)
    else:
        segments = parse_option(options.segments)

    M = parse_option(options.num_components)

    dict_nodes_rem = {}
    node = 0
    for sh in shanks:
        for sg in segments:
            for nc in M:   
                for f in range(0, int(options.num_folds)+1):   
                    for rep in range(1, int(options.num_replicas)+1):
                        # NOTE: f == 0 is interpreted as num_folds == 1
                        node += 1
                        dict_nodes_rem[node] = [int(sh), int(sg),\
                                                   int(nc), int(f),\
                                                   int(rep)]
    return (node, dict_nodes_rem)

# ----------------------------------------------------------------------

def map_nodes_extract_best_solution(ft_file_path, options):
    if options.shanks == "all":
        num_shanks = get_num_shanks_from_data(ft_file_path)
        shanks = range(1, num_shanks+1)
    else:
        shanks   = parse_option(options.shanks)
    if options.segments == "all":
        num_segments = get_num_segments_from_data(ft_file_path)
        segments = range(1, num_segments+1)
    else:
        segments = parse_option(options.segments)
        
    M  = parse_option(options.num_components)
    nf =          int(options.num_folds)

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

def map_nodes_cross_validation(ft_file_path, options):
    if options.shanks == "all":
        num_shanks = get_num_shanks_from_data(ft_file_path)
        shanks = range(1, num_shanks+1)
    else:
        shanks   = parse_option(options.shanks)
    if options.segments == "all":
        num_segments = get_num_segments_from_data(ft_file_path)
        segments = range(1, num_segments+1)
    else:
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

def map_nodes_best_num_components(ft_file_path, options):
    if options.shanks == "all":
        num_shanks = get_num_shanks_from_data(ft_file_path)
        shanks = range(1, num_shanks+1)
    else:
        shanks   = parse_option(options.shanks)
    if options.segments == "all":
        num_segments = get_num_segments_from_data(ft_file_path)
        segments = range(1, num_segments+1)
    else:
        segments = parse_option(options.segments)

    dict_nodes_best_num_components = {}
    node = 0
    for sh in shanks:
        for sg in segments:
            node += 1
            dict_nodes_best_num_components[node] = [int(sh), int(sg)]
    return (node, dict_nodes_best_num_components)

# ----------------------------------------------------------------------

def get_num_segments_from_data(ft_file):
    mat = scipy.io.loadmat(ft_file,squeeze_me=True,struct_as_record=False)
    try:
        num_segments = len(mat['featStruct'].__dict__['seg'])
    except:
        sys.exit("Please, explicitly specify the segments to be processed (with option -T)")
    return num_segments

# ----------------------------------------------------------------------

def get_num_shanks_from_data(ft_file):
    mat = scipy.io.loadmat(ft_file,squeeze_me=True,struct_as_record=False)
    try:
        num_shanks = len(mat['featStruct'].__dict__['seg'][0].__dict__['shank'])
    except:
        sys.exit("Please, explicitly specify the shanks to be processed (with option -S)")
    return num_shanks

# ----------------------------------------------------------------------
if __name__ == "__main__":
   
    usage = "Usage: \n\
    %prog input_*ft.mat_file(s) [options (-h to list)]"   
    # method must be one of: 'ms','em'
    # node # will be mapped to [segment#, shank#, #components]

    parser = optparse.OptionParser(usage=usage, version="%%prog ")
    parser = spike_sorting_command_line_parser(parser)
    (options, args) = parser.parse_args()

    if options.compile or options.compile_all:
        compile_code(options)
#   print "len(args)=", len(args)
    if len(args) > 0:
        input_files = []
        for arg in args:
            input_files1 = parse_input_data(arg, tonic_data)
            for if1 in input_files1:
                if not if1 in input_files:
                    input_files.append(if1)
        print "input files=", input_files
        if options.method == "kk":
            options.num_replicas = 1;
        for input_file in input_files:
            ft_file_path = os.path.join(tonic_data, input_file)
            num_nodes_rem, dict_nodes_rem = map_nodes_rem(ft_file_path, options)
            if options.segments == "all":
                num_segments = get_num_segments_from_data(input_file)
            else:
                num_segments = len(parse_option(options.segments))
            if options.shanks == "all":
                num_shanks = get_num_shanks_from_data(input_file)
            else:
                num_shanks = len(parse_option(options.shanks))
            if options.verbose:
                print "num_segments=", num_segments, " num_shanks=", num_shanks    
            print "len(args)=", len(args), " args=", args 
            if len(options.inode) == 0:
                sort_spikes_high_level(ft_file_path, options.method, num_nodes_rem, options)
                if options.verbose:
                    print "num_nodes_rem= ", len(dict_nodes_rem.keys())
            else:
                print "options.inode=", options.inode 
                node = int(options.inode)
                rem_low_level(ft_file_path, node, dict_nodes_rem, options)
    elif options.compile or options.compile_all:
        sys.exit(2)
    else:
        parser.print_usage()
        sys.exit(2)
