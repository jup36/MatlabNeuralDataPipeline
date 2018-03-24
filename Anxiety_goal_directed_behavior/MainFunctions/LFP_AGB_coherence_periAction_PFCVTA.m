%This script uses the chronux function 'cohgramc' to get the coherence
% between PFC and VTA LFP timeseries in the peri-event time windows.   

clear all; clear functions; clc

%%
addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab\Functions'));

fileDirectory = 'C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT\LFP'; cd(fileDirectory);   % directory for raw file
load('LFP_AGB_PSD_PFC_092815','raw','tstpbeh','lfpID');   % load PFC raw lfp time series 
PFCraw = raw; PFClfpID = lfpID; clearvars raw lfpID;

load('LFP_AGB_PSD_VTA_092815','raw','lfpID');   % load VTA raw lfp time series 
VTAraw = raw; VTAlfpID = lfpID; clearvars raw lfpID;

params.Fs=1000; % sampling frequency
params.fpass=[0 100]; % band of frequencies to be kept
params.tapers=[5 9]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[2 0.05];
params.trialave=1;
movingwin=[0.5 0.05];
% segave=1;
% direction=5;

% wintrig=[5*movingwin(1) 5*movingwin(1)];
% winseg=2*movingwin(1);

% preallocate cells
C = cell(size(VTAraw.ispk,3),1); % cell containing each animal's trial-to-trial coherence for event time series
Cbase = cell(size(VTAraw.ispk,3),1); % cell containing each animal's trial-to-trial coherence for base time series
rawC.B1 = cell(size(VTAraw.ispk,3),1); rawC.B2 = cell(size(VTAraw.ispk,3),1); rawC.B3 = cell(size(VTAraw.ispk,3),1); % cell containing each animal's raw coherence averaged across trials in each block 
normC.B1 = cell(size(VTAraw.ispk,3),1); normC.B2 = cell(size(VTAraw.ispk,3),1); normC.B3 = cell(size(VTAraw.ispk,3),1); % cell containing each animal's normalized coherence averaged across trials in each block

%% Compute raw and normalized coherence data
for f=1:size(VTAraw.ispk,3) % # of files (animals)
    
    PFClfpidx = zeros(length(PFClfpID),1); % index for the matching PFC file
    % find the matching PFC lfp data for the current VTA lfp data
    for o = 1:length(PFClfpID) % the number of PFClfp files; As there are more PFC files recorded, get VTA lfp first, and find PFC lfp that matches to that of VTA!
        PFClfpidx(o,1) = strcmp(VTAlfpID{f}(1:11),PFClfpID{o}(1:11));    % identify the valid tstpbeh by matching names
    end
    clearvars o
    PFClfpidx = logical(PFClfpidx);
    
    valtridx = zeros(size(tstpbeh,2),1); % index for the matching tstpbeh, required to get the currently valid trial number
    for o = 1:size(tstpbeh,2) % the number of PFClfp files; As there are more PFC files recorded, get VTA lfp first, and find PFC lfp that matches to that of VTA!
        valtridx(o,1) = strcmp(VTAlfpID{f}(1:9),tstpbeh(o).name(1:9));    % identify the valid tstpbeh by matching names
    end
    clearvars o
    currnumbtr = length(tstpbeh(find(valtridx==1)).ispk);    % get the number of the current valid trials
      
    % count valid trials across all trials
    valtr = nan(size(VTAraw.ispk(:,1,f),1),1); % to count valid trials
    
    for t = 1:currnumbtr    % # of trials
        if f==1 && t==1     % using the very 1st trial of the 1st animal
           [~,~,~,~,~,time,freq,~,~,~]=cohgramc(PFCraw.ispk{t,1,PFClfpidx}',VTAraw.ispk{t,1,f}',movingwin,params); % get the # of time and freq bins 
           [~,~,~,~,~,basetime,basefreq,~,~,~]=cohgramc(PFCraw.base{t,1,PFClfpidx}',VTAraw.base{t,1,f}',movingwin,params); % get the # of time and freq bins  
        end
        
        % Multi-taper time-frequency coherence for event timeseries (cohgramc)
        if isnan(VTAraw.ispk{t,1,f})|isnan(PFCraw.ispk{t,1,PFClfpidx})==1;  % invalid trials 
            valtr(t) = 0;
            C{f}(:,:,t)=nan(length(time),length(freq));   % put NaN for invalid trials
        else    % valid trials
            valtr(t) = 1;
            [C{f}(:,:,t)]=cohgramc(PFCraw.ispk{t,1,PFClfpidx}',VTAraw.ispk{t,1,f}',movingwin,params);   % get the trial-to-trial coherence between PFC and VTA LFP segments 
        end
        
        % Multi-taper time-frequency coherence for base timeseries (cohgramc)
        if isnan(VTAraw.base{t,1,f})|isnan(PFCraw.base{t,1,PFClfpidx})==1;  % invalid trials 
            Cbase{f}(:,:,t)=nan(length(basetime),length(basefreq));   % put NaN for invalid trials
        else
            [Cbase{f}(:,:,t)]=cohgramc(PFCraw.base{t,1,PFClfpidx}',VTAraw.base{t,1,f}',movingwin,params);   % get the trial-to-trial coherence between PFC and VTA LFP segments 
        end       
        fprintf('processed trial #%d\n', t);    % report the progress (completion of each valid trial)
    end
    clearvars t
       
    if nansum(valtr) > 150/2;  % if more than half trials are valid, include the coherence data for rawC and normC    
       rawC.B1{f} = C{f}(:,:,1:50); rawC.B2{f} = C{f}(:,:,51:100); rawC.B3{f} = C{f}(:,:,101:end); 
       
       tempavgbasecoh = nanmean(nanmean(Cbase{f},1),3);     % mean baseline coherence (averaged across time bins and trials)
       tempstdbasecoh = nanstd(nanstd(Cbase{f},0,1),0,3);   % std baseline coherence (averaged across time bins and trials)
       
       tempavgbasecohrep = repmat(tempavgbasecoh,length(time),1,50);
       tempstdbasecohrep = repmat(tempstdbasecoh,length(time),1,50);
       
       normC.B1{f} = (C{f}(:,:,1:50)-tempavgbasecohrep)./tempstdbasecohrep;     % normalization to the baseline (total trials) coherence
       normC.B2{f} = (C{f}(:,:,51:100)-tempavgbasecohrep)./tempstdbasecohrep;   % normalization to the baseline (total trials) coherence
       normC.B3{f} = (C{f}(:,:,101:end)-tempavgbasecohrep(:,:,1:size(C{f}(:,:,101:end),3)))./tempstdbasecohrep(:,:,1:size(C{f}(:,:,101:end),3));  % normalization to the baseline (total trials) coherence
       
       normC.B1avg(:,:,f) = nanmean(normC.B1{f},3); % average across B1 trials
       normC.B2avg(:,:,f) = nanmean(normC.B2{f},3); % average across B2 trials
       normC.B3avg(:,:,f) = nanmean(normC.B3{f},3); % average across B3 trials
       
       rawC.B1avg(:,:,f) = nanmean(rawC.B1{f},3);   % average across B1 trials
       rawC.B2avg(:,:,f) = nanmean(rawC.B2{f},3);   % average across B2 trials
       rawC.B3avg(:,:,f) = nanmean(rawC.B3{f},3);   % average across B3 trials
       
    else
        
    end  
    fprintf('processed file #%d\n', f);    % report the progress (completion of each VTA-PFC pair)
end

save('PFC_VTA_periAction_coherence','C','Cbase','normC','rawC','time','freq') % save the raw and normalized coherence data

%% plot the normalized coherence averaged across animals
B1avgnormC = nanmean(normC.B1avg,3);  % Block1 Coherence average across animals
B2avgnormC = nanmean(normC.B2avg,3);  % Block2 Coherence average across animals
B3avgnormC = nanmean(normC.B3avg,3);  % Block3 Coherence average across animals

figure; % Block1 Baseline-normalized average coherence
imagesc(imrotate(B1avgnormC,90)) % imrotate at 90 degrees counter-clockwise
xlim([20 60])   % -1 to 1 sec
set(gca,'xtick',30:10:50); % xTick -0.5 to 0.5
set(gca,'ytick',0:40:200); % yTick 0 to 100 Hz with steps of 20 Hz
set(gca,'TickDir','out');
set(gca,'FontSize',14);
set(gca,'Layer','top');
set(gca,'GridLineStyle','none');
colorbar('EastOutside');
caxis([-1 2])

cd('C:\Documents and Settings\jup36\Desktop\Matlab\colorMaps');
load('hotcolormap');
colormap(hotcolormap);
clear custom*

xlabel('Time (sec)','FontSize',14)
ylabel('Frequency (Hz)','FontSize',14)   

figure; % Block2 Baseline-normalized average coherence
imagesc(imrotate(B2avgnormC,90)) % imrotate at 90 degrees counter-clockwise
xlim([20 60])   % -1 to 1 sec
set(gca,'xtick',30:10:50); % xTick -0.5 to 0.5
set(gca,'ytick',0:40:200); % yTick 0 to 100 Hz with steps of 20 Hz
set(gca,'TickDir','out');
set(gca,'FontSize',14);
set(gca,'Layer','top');
set(gca,'GridLineStyle','none');
colorbar('EastOutside');
caxis([-1 2])

cd('C:\Documents and Settings\jup36\Desktop\Matlab\colorMaps');
load('hotcolormap');
colormap(hotcolormap);
clear custom*

xlabel('Time (sec)','FontSize',14)
ylabel('Frequency (Hz)','FontSize',14)   

figure; % Block3 Baseline-normalized average coherence
imagesc(imrotate(B3avgnormC,90)) % imrotate at 90 degrees counter-clockwise
xlim([20 60])   % -1 to 1 sec
set(gca,'xtick',30:10:50); % xTick -0.5 to 0.5
set(gca,'ytick',0:40:200); % yTick 0 to 100 Hz with steps of 20 Hz
set(gca,'TickDir','out');
set(gca,'FontSize',14);
set(gca,'Layer','top');
set(gca,'GridLineStyle','none');
colorbar('EastOutside');
caxis([-1 2])

cd('C:\Documents and Settings\jup36\Desktop\Matlab\colorMaps');
load('hotcolormap');
colormap(hotcolormap);
clear custom*

xlabel('Time (sec)','FontSize',14)
ylabel('Frequency (Hz)','FontSize',14)   



