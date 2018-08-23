clear all;
close all;
clc;
cd('E:\NEX\ADQ\FR0108')
% cd('S:\MoghaddamLab\Wood\Meredyth Data\FR0101')
% cd('S:\MoghaddamLab\Wood\Meredyth Data\Session1')
addpath(genpath('C:\Program Files\MATLAB\chronux'))
%% user defined variables
uv.reg1 = [1:8];
uv.reg2 = [9:16];
uv.dT = 0.025;                                                              %time bin size in sec
uv.selEvent = 'EVT13';                                                      %event of interest
uv.evtWin = [10 10];                                                        %time range before and after event in sec
uv.baseEvent = 'EVT10';                                                     %event that defines baseline (end/start of ITI)
uv.baseWin = [8 -4];                                                        %baseline period relative to ITI event
uv.time = [-1*uv.evtWin(1,1):uv.dT:uv.evtWin(1,2)];                         %time scale for neuronal activity
uv.smoothWidth = 3;                                                         %number of bins to smooth data by
%% get info about files in directory
files = dir('MW*');

%% load spike and event data
for f = 1:size(files,1)                                                     %for all files
    fName = files(f,1).name;                                                
    %% get info about each file
    [~,names,types]  =  nex_info(fName);                                    %get the file's info                                   
    evts = types(:,1) == 1;                                                 %logical for event type data
    units = types(:,1) == 0;                                                %logical for unit timestamp type data
    wforms = types(:,1) == 3;                                               %logical for waveform type data
    unitNames = names(units,:);                                             %get the names of each unit
    evtNames = names(evts,:);                                               %get the names of each event
    wfNames = names(wforms,:);                                              %get the name of each waveform
    chans{f,1} = str2num(unitNames(:,5:6));                                 %convert unit name to numerical val
    clear types units evts names wforms
    
    %% get all spike times per session
    for i = 1:size(unitNames,1)                                             %for every unit in the current file
        [~,st{f,i}]= nex_ts(fName,unitNames(i,:));                          %load in the neuronal time stamps 
    end
    
    %% get the event timestamps
    [~,evt{f,1}]= nex_ts(fName,uv.baseEvent);                               %load in the ITI edge timestamps
    [~,evt{f,2}]= nex_ts(fName,uv.selEvent);                                %load in the event of interest timestamps
    clear unitNames evtNames i idx
    %% get waveforms
    for w = 1:size(wfNames,1)
        [adfreq,~,~,~,waveforms{f,w}]= nex_wf(fName,wfNames(w,:));
    end
end
chanNum = cell2mat(chans);                                                  %convert cell format to matrix
regIDX = chanNum <= uv.reg1(1,end);                                         %logical index for channels on first array
clear f chanNum w

%% align spike times relative to event times
counter = 0;                                                            
for sess = 1:size(evt,1)                                                    %for each session
    [maxU,~] = cellfun(@size,chans(sess,:),'UniformOutput',1);              %find how many units are in the current session
    for u = 1:maxU
        counter = counter + 1;                                              %increment counter (used to store data by rows)
        spikes = cell2mat(st(sess,u));                                      %convert to matrix format the spike times
        base = cell2mat(evt(sess,1));                                       %''          ''               base edge times
        events = cell2mat(evt(sess,2));                                     %''          ''               event times
        
        %% get spike times aligned to event
        trigST(1,:) = createdatamatpt(spikes,base,uv.baseWin);              %align spike times to the beginning of each event window in each trial
        trigST(2,:) = createdatamatpt(spikes,events,uv.evtWin);                     
        clear base spikes events
        
        %%
        for ev = 1:size(trigST)                                             %for each event (baseline and event of interest)
            %% initialize matrix with NaNs
            tempTrialMat = NaN(size(trigST,2),1);                           %initialize a matrix of NaN values
            
            %% get spike times in analysis window and align to event
            for trial = 1:size(trigST(ev,:),2)                              %for each trial
                if length(trigST(ev,trial).times) > size(tempTrialMat,2)    %if the length of the trial matrix is less than the length of the spike time vector
                    lenDiff = diff([size(tempTrialMat,2),...                %get the difference between the number of cols in the matrix and the len of the vector
                        length(trigST(ev,trial).times)]);
                    tempTrialMat(:,end+1:end+lenDiff) = NaN;                %append the matrix with NaNs
                end
                tempTrialMat(trial,1:length(trigST(ev,trial).times))=...
                    trigST(ev,trial).times - uv.evtWin(1,1);                %add the spike times to the matrix
            end                                                                                                           
            clear spikes events trial lenDiff 
            
            %% get spike counts and firing rates
            count = histc(tempTrialMat,uv.time,2);                          %spike counts binned at dT (same format as above)
            count = count(:,1:end-1);                                       %remove the last bin (which is always empty)
            trialsFR{counter,ev} = count / uv.dT;                           %divide the spike count by binsize to yield firing rate
            FR(ev,1).raw(counter,:) = mean(count / uv.dT)
            clear count tempTrialMat
        end
        clear ev
    end
    clear u maxU trigST
end
clear sess counter chans 

%% Z score norm the data
for u = 1:size(trialsFR,1)                                                  %for every unit
    %% get mean and SD
    base = cell2mat(trialsFR(u,1));                                         %get the baseline trials
    sdBaseFR(u) = std(base(:));                                             %get the whole window-wide standard deviation firing rate
    avBaseFR(u) = mean(base(:));                                            %get the whole window-wide mean firing rate
    
    %% z score, smooth, and place into matrix
    for ev = 1:size(evt,2)                                                  %for all events
        temp = bsxfun(@minus,cell2mat(trialsFR(u,ev)),avBaseFR(u));         %subtract off the average
        normTrials{u,ev} = temp/sdBaseFR(u);                                %divide by sd (Z = (x - mean)/sd)
        smNormTrials{u,ev} = smoothts(cell2mat(normTrials(u,ev)),'b',uv.smoothWidth); %smooth
        FR(ev,1).norm(u,:) = mean(normTrials{u,ev});                        %get the normalized mean firing rate across all trials
        FR(ev,1).smNorm(u,:) = mean(smNormTrials{u,ev});                    %get the smoothed and normalized mean firing rate across all trials
    end
    clear base ev temp
    %%
end
clear u

%% divide data by region
for ev = 1:size(FR,1)                                                       %for all events
    reg1FR(ev,1).norm = FR(ev,1).norm(regIDX,:)                             %reg 1 data 
    reg1FR(ev,1).smNorm = FR(ev,1).smNorm(regIDX,:)
    reg2FR(ev,1).norm = FR(ev,1).norm(~regIDX,:)                            %reg 2
    reg2FR(ev,1).smNorm = FR(ev,1).smNorm(~regIDX,:)
end

%% get baseline firing rate
BaseFR=FR(1,1);
BaseFR = BaseFR.raw;
BaseFR = BaseFR(:,1:200);
avbFR= mean (BaseFR,2);
%% manually select the data to plot

data = reg1FR(2,1).smNorm;                                                  %select the data to plot
av = mean(data);                                                            %average across units
se = std(data)/sqrt(size(data,1));                                          %compute the se across units
t  = uv.time(1:end-1);

fill([t,fliplr(t)],[av+se,fliplr(av-se)],'r')                               %fill plot
xlim([-.5 .5])
xlabel('Time (sec)')
ylim([-0.2 4])
ylabel('Normalized Population Activity (Z Score)')
hold on
plot(t,av)
hold on

data_2 = reg2FR(2,1).smNorm;                                                  %select the data to plot
av_2 = mean(data_2);                                                            %average across units
se_2 = std(data_2)/sqrt(size(data_2,1));                                          %compute the se across units
t_2  = uv.time(1:end-1);                                              

fill([t_2,fliplr(t_2)],[av_2+se_2,fliplr(av_2-se_2)],'b')                               %fill plot
xlim([-.5 .5])
xlabel('Time (sec)')
ylim([-0.2 4])
ylabel('Normalized Population Activity (Z Score)')
hold on
plot(t_2,av_2)


%%epoch isolation, consolidation, p-value

V = data(:,400:409);                                                        %isolate ten time bins (0.25s)around the event
S = data_2(:,400:409);

VTA = sum(V,2);                                                             %number of units in each region
SNc = sum(S,2);
% % 
% % for u = 1:size(S,1);                                                        %%Find peak firing rate for each unit in SNc
% %     temp = max(S(u,:));
% %     peakS(u,1) = temp;
% %     clear temp
% % end
% % 
% % for v = 1:size(V,1);                                                        %%Find peak firing rate for each unit in VTA
% %     temp = max(V(v,:));
% %     peakV(v,1) = temp;
% %     clear temp
% % end
% % 
% % 
% % [h,p] = ttest2(peakV,peakS)                                                 %%unpaired ttest for peak firing rate of VTA and SNc from the middle ten bins of that event window

[h,p] = ttest2(VTA,SNc)

