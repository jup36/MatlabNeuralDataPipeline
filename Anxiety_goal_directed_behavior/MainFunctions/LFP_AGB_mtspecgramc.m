%This script is to get the raw and normalized power spectral density (psd)
% values from the continuous LFP recording data. To run this script, NEX
% file, tstp (timestamps of events) must be ready for each animal. Also,
% the valid LFP channel of each animal must be specified in the channelList to run nex_cont. 
% 'eventpsdAGB' is the function to run 'mtspecgramc' for each trial. 
% 'mtspecgramc' is used for FFT with moving-window and multitaper. 
% 'createdatamatc' is used to get the peri-event time windows. 
% JP11 has been excluded from analysis

clear all; clear functions; clc

%%
addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab\Functions'));

%% load signal and timestamps modify this!
filedirectory = 'C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\NEX';
savedirectory = 'C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\NEX';
cd(filedirectory);
load('SUA_AGB_TSTP_120513');                    % load timestamps, the timestamps are in the unit of sec                  
fileList = dir('AGB_LFP_VTA*');
channelList = ['FP20';'FP18';'FP18'];           % channel ID, JP12, 14, 15 (JP11 has been excluded from analysis) 
saveName = 'AGB_VTA_psd_120513.mat';

%% set parameters
sv.timewindow = 500;        % msec
sv.step = 100;              % msec
sv.baseendpoint = -0.5;     % sec, baseline period will start -3.5 to -0.5 sec from the onset of each cue           
params.trialave = 0;
params.pad = -1;
movingwin = [0.5 0.1];      % msec
params.fpass = [0 100];     % 0 to 100 Hz
params.tapers = [3 5];      % the number of tapes (5 is default)
params.Fs = 1000;           % sampling rate (1000 Hz) 

%% get psd(0-100Hz) around the whole responses (without classifying the trials) in each session of each animal  
for i = 1:size(fileList,1)        % number of files or animals in this dose group
    
    % get the animal ID, to get the correct timestamps for the current data
    switch fileList(i).name(1,13:16)
        case 'JP11'
            animalID = 1;
        case 'JP12'
            animalID = 2;
        case 'JP14'
            animalID = 3;
        case 'JP15'
            animalID = 4;
    end
    
    cue_ts = tstp.cue{cell2mat(tstp.cue(:,2)) == animalID,1};                 % get the current cue timestamps by matching the animal ID
    ispk_ts = tstp.ispk{cell2mat(tstp.ispk(:,2)) == animalID,1};              % get the current ispk timestamps by matching the animal ID
    pellet_ts = tstp.pellet{cell2mat(tstp.pellet(:,2)) == animalID,1};        % get the current pellet timestamps by matching the animal ID
    iti_ts = tstp.iti{cell2mat(tstp.iti(:,2)) == animalID,1};                 % get the current iti timestamps by matching the animal ID
    
    [~,~,~,~,d] = nex_cont(fileList(i,1).name,channelList(i,:));              % get the continuous data of each animal 
    
    meand = mean(d,2);      % get the mean to set the cut-off criteria
    stdd = std(d,0,2);      % get the std to set the cut-off criteria
    % plot(d)                    % use this to plot the raw LFP voltage signal 
    % axis([0 size(d,2) -1 1])   % scaling and appearance of the plot    
    %% this part is to get the upper and lower limits to use as cut-off of the raw signal 
    upperlimit = meand + 2*stdd;
    lowerlimit = meand - 2*stdd;
        
    tempbase = cell(1,1);   % tempbase contains valid baselines  
    tempbasecount = 0;      % counter for valid tempbase
          
    %% loop to get peri-evt (ISPK) psd using the function longeventpsd containing the chronux function mtspectramc    
    for ii = 1:size(ispk_ts,1)-1          % number of total ISPKs - 1 (# of trials), because the last trial was most likely to be truncated 
        [ psd.ispk{ii,1,i}, raw.ispk{ii,1,i}, psd.norm_ispk{ii,1,i}, psd.base_ispk{ii,1,i}, raw.base_ispk{ii,1,i}, t, f ] = eventpsdAGB( d, ispk_ts(ii,1), cue_ts, movingwin, params, sv, upperlimit, lowerlimit );       % get t and f to use it for later in the script
        if sum(isnan(psd.base_ispk{ii,1,i}))==0
           tempbasecount = tempbasecount + 1;
           tempbase{tempbasecount,1}(:,:) = psd.base_ispk{ii,1,i};
        else 
        end
     
    end
    clearvars ii j
    
    % normalization to the avg collected baseline
    [psd.normavg_ispk(1:size(ispk_ts,1)-1,1,i), collectbasemat(i,:), collectbasestdmat(i,:)] = psdnormalization(psd.ispk(1:size(ispk_ts,1)-1,:,i), tempbase, t, f);       % number of total ISPKs - 1 (# of trials), because the last trial was most likely to be truncated 
    
    %% loop to get peri-evt (CUE) psd using the function longeventpsd containing the chronux function mtspectramc    
    for ii = 1:size(cue_ts,1)-1         % number of total cues - 1 (# of trials), because the last trial was most likely to be truncated 
        [ psd.cue{ii,1,i}, raw.cue{ii,1,i}, psd.norm_cue{ii,1,i}, psd.base_cue{ii,1,i}, raw.base_cue{ii,1,i}, ~, ~ ] = eventpsdAGB( d, cue_ts(ii,1), cue_ts, movingwin, params, sv, upperlimit, lowerlimit );       % get t and f to use it for later in the script
        if sum(isnan(psd.base_cue{ii,1,i}))==0
           tempbasecount = tempbasecount + 1;
           tempbase{tempbasecount,1}(:,:) = psd.base_cue{ii,1,i};
        else 
        end
     
    end
    clearvars ii j
    
    % normalization to the avg collected baseline
    [psd.normavg_cue(1:size(cue_ts,1)-1,1,i), collectbasemat(i,:), collectbasestdmat(i,:)] = psdnormalization(psd.cue(1:size(cue_ts,1)-1,:,i), tempbase, t, f);       % number of total cues - 1 (# of trials), because the last trial was most likely to be truncated %% loop to get peri-evt (CUE) psd using the function 

    %% loop to get peri-evt (PELLET) psd using the function longeventpsd containing the chronux function mtspectramc    
    for ii = 1:size(pellet_ts,1)-1          % number of total pellets - 1 (# of trials), because the last trial was most likely to be truncated 
        [ psd.pellet{ii,1,i}, raw.pellet{ii,1,i}, psd.norm_pellet{ii,1,i}, psd.base_pellet{ii,1,i}, raw.base_pellet{ii,1,i}, ~, ~ ] = eventpsdAGB( d, pellet_ts(ii,1), cue_ts, movingwin, params, sv, upperlimit, lowerlimit );       % get t and f to use it for later in the script
        if sum(isnan(psd.base_pellet{ii,1,i}))==0
           tempbasecount = tempbasecount + 1;
           tempbase{tempbasecount,1}(:,:) = psd.base_pellet{ii,1,i};
        else 
        end
     
    end
    clearvars ii j
    
    % normalization to the avg collected baseline
    [psd.normavg_pellet(1:size(pellet_ts,1)-1,1,i), collectbasemat(i,:), collectbasestdmat(i,:)] = psdnormalization(psd.pellet(1:size(pellet_ts,1)-1,:,i), tempbase, t, f);       % number of total pellets - 1 (# of trials), because the last trial was most likely to be truncated 
    
    %% loop to get peri-evt (ITI) psd using the function longeventpsd containing the chronux function mtspectramc    
    for ii = 1:size(iti_ts,1)-1         % number of total itis - 1 (# of trials), because the last trial was most likely to be truncated 
        [ psd.iti{ii,1,i}, raw.iti{ii,1,i}, psd.norm_iti{ii,1,i}, psd.base_iti{ii,1,i}, raw.base_iti{ii,1,i}, ~, ~ ] = eventpsdAGB( d, iti_ts(ii,1), cue_ts, movingwin, params, sv, upperlimit, lowerlimit );       % get t and f to use it for later in the script
        if sum(isnan(psd.base_iti{ii,1,i}))==0
           tempbasecount = tempbasecount + 1;
           tempbase{tempbasecount,1}(:,:) = psd.base_iti{ii,1,i};
        else 
        end
     
    end
    clearvars ii j
    
    % normalization to the avg collected baseline
    [psd.normavg_iti(1:size(iti_ts,1)-1,1,i), collectbasemat(i,:), collectbasestdmat(i,:)] = psdnormalization(psd.iti(1:size(iti_ts,1)-1,:,i), tempbase, t, f);       % number of total itis - 1 (# of trials), because the last trial was most likely to be truncated 
    
end
clearvars i meand stdd upperlimit lowerlimit d

%%
cd(savedirectory);
clearvars filedirectory savedirectory;
save(saveName);
