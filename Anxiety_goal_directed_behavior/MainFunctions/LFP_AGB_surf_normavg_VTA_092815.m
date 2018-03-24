%This script is to collect and display the normavg psd values of the entire
% correct trials in all valid animals. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input file: SUA_AGB_tstpbeh, LFP_AGB_PSD_VTA                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clear functions; clc

%%
addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab\Functions'));

%% load signal and timestamps modify this!
fileDirectory1 = 'C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT\LFP'; cd(fileDirectory1);   % directory for psd file
filename = 'LFP_AGB_PSD_VTA_092815';    % the file containing event-related PSDs  

load(filename,'psd','raw','t','f','params','sv','tstpbeh','lfpID');   
freq = f;

%% user variables, modify this!!
uv.block = 3;

%% collect valid trials of all the animals  
valsufftrcount = 1;     % valid and sufficient number of trials 
valsuffanimal = 1;

eachvalpsdcell = cell(size(psd.ispk,3),1); % cell to contain the psd data
for f = 1:size(psd.ispk,3)   % # of animals (files)
    
    % count valid trials across all trials
    valtr = nan(size(psd.normavg_ispk(:,1,f),1),1);
    for t = 1:size(psd.normavg_ispk(:,1,f),1)
        if isnan(psd.normavg_ispk{t,1,f})==1;
            valtr(t) = 0;
        else
            valtr(t) = 1;
        end
    end
    
    switch uv.block
        case 1
           if lfpID{f}(11)=='A'
               crtrials = 1:50;
           elseif lfpID{f}(11)=='D'
               crtrials = 101:150; end
        case 2
           crtrials = 51:100; 
        case 3  
           if lfpID{f}(11)=='A'
               crtrials = 101:150;
           elseif lfpID{f}(11)=='D'
               crtrials = 1:50; end
    end
    
    if sum(valtr) >= round(size(psd.normavg_ispk(:,1,f),1)/2);  % if more than half trials are valid
        temppsd = psd.normavg_ispk(crtrials,1,f);       % index the psd.normavg_evt of trials of interest
        valtrcount = 0; % reset the # of valid trials for new animals
        % loop to collect valid trials of each animal
        for t = 1:size(temppsd,1)      % # of trials of interest
            if isnan(temppsd{t,1}) == 1
            else
                if isempty(temppsd{t,1}) == 1      % empty trial: nothing happens
                else    % a valid trial: proceed!
                    valtrcount = valtrcount + 1;    % count the valid trial
                    tempvalpsd(:,:,valtrcount) = temppsd{t,1};    % put the psd mat of the valid trial to the temppsd, (9:29,1:31) corresponds to -1 to 1 sec, 0 to 60 Hz
                end
            end
        end
        clearvars t
        
        if valtrcount < 8      % in case, there is no valid trials at all
            eachvalpsdcell{f,1} = NaN;
            eachnumbvaltr(f,1) = valtrcount;
        elseif valtrcount >= 8  % in case, there is 8 or more valid trials
            eachavgvalpsd(:,:,valsuffanimal) = nanmean(tempvalpsd,3);        % get the average across trials for EACH ANIMAL and put the data to the eachavgvalpsd
            eachvalpsdcell{f,1} = tempvalpsd;
            eachnumbvaltr(f,1) = valtrcount;
            eachvalsuffmat(:,:,valsufftrcount:valsufftrcount+valtrcount-1) = tempvalpsd;         % eachvalsuffmat collects all the valid trials of interest, of which number of trials exceeding 8
            valsufftrcount = valsufftrcount + valtrcount;
            valsuffanimal =  valsuffanimal + 1;
        end
        clearvars tempvalpsd
        
    else % if no more than half trials are valid
         eachvalpsdcell{f,1} = NaN; % if no more than half trials are valid
    end
        
end
clearvars f temp*

%% select files to be plotted
plotvalpsdmat = nanmean(eachavgvalpsd(:,:,1:end),3);    % average across rats 
%plotvalpsdmat = nanmean(eachvalpsdcell{25},3);

%% smoothing for the averaged psd in both time and frequency dimensions
for i = 1:size(plotvalpsdmat,2);          % for smoothing in time domain
    sm_avgvalpsdmat(:,i) = smooth(plotvalpsdmat(:,i),5);      % smoothing across row (time)
end
clearvars i

%% image using surf figure
[n_row, n_col]=size(sm_avgvalpsdmat);
time_bin=4;      % 4 sec, the length of time window     
max_freq=100;

% time_scale = (0+time_bin/n_row)-2:time_bin/n_row:(time_bin)-2;             % the length of the entire timewindow is 4 sec, which is divided into 2 + 2 around the event. 
% 
% [x,y]=meshgrid(time_scale,[max_freq/n_col:max_freq/n_col:max_freq]);
% 
% x = x';         % to match dimension
% y = y';         % to match dimension

figure;
%imagesc(sm_avgvalpsdmat)
imagesc(imrotate(sm_avgvalpsdmat(:,3:51),90))
% Rotate First!!
xlim([20 60])
set(gca,'xtick',30:10:50);
set(gca,'TickDir','out');
set(gca,'FontSize',14);
set(gca,'Layer','top');
set(gca,'GridLineStyle','none');
colorbar('EastOutside');
caxis([-1.5 1.5])
% figure = surf(x,y,sm_avgvalpsdmat);
% set(figure,'edgecolor','none')
% view(0,90);
cd('C:\Documents and Settings\jup36\Desktop\Matlab\colorMaps');
load('hotcolormap');
colormap(hotcolormap);
clear custom*

xlabel('Time (sec)','FontSize',14)
ylabel('Frequency (Hz)','FontSize',14)   

% save figure
% cd(savefigdirectory);
% figurename = [figurename];
% print('-depsc','-tiff','-r500',figurename)