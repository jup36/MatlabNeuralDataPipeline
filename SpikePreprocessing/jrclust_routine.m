%SPIKE SORTING ROUTINE USING JRC
% Please make sure the data folder has '*.ap.bin' file, and proper '*.prb'
% file. 

clear all; clear functions; clc

addpath(genpath('/Volumes/RAID2/parkj/MATLAB/JRCLUST'));  % jrclust folder
cd('/Volumes/RAID2/parkj/NeuralData/IT03_Ldms_M1_013118/Matfiles') % data folder
%jrc git-pull % run for jrc updates

jrc makeprm ITPhys03_180131_S023_m1dstr_stim_m1_g1_t0.imec.ap.bin imec3_opt3.prb % create the parameter file

jrc spikesort ITPhys03_180131_S023_m1dstr_stim_m1_g1_t0.imec.ap_imec3_opt3.prm   % cluster automatically 

jrc manual ITPhys03_180131_S023_m1dstr_stim_m1_g1_t0.imec.ap_imec3_opt3.prm

%jrc3 probe imec3_opt3.prb  % show the probe layout

%jrc3 traces ITPhys03_180131_S023_m1dstr_stim_m1_g1_t0.imec.ap_imec3_opt3.prm  % plot traces

%jrc3 detect ITPhys03_180131_S023_m1dstr_stim_m1_g1_t0.imec.ap_imec3_opt3.prm  % detect spikes

%jrc3 describe ITPhys03_180131_S023_m1dstr_stim_m1_g1_t0.imec.ap_imec3_opt3.prm

% export csv - no need to export csv files anymore (since 2/22/18)! Subsequent analyses
% will just use information from '_jrc.mat'.  
%cd('/Volumes/RAID2/parkj/NeuralData/WR26_Lacc_dms_101917/Matfiles') % RAID directory mac

%jrc3 exportcsv ITPhys03_180131_S023_m1dstr_stim_m1_g1_t0.imec.ap_imec3_opt3.prm    % spike timing, cluster number and max. site locations are saved in a tabular format
