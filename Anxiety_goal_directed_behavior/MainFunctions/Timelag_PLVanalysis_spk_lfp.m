%This script is to perform the time-lag Phase-locking analysis between VTA 
% spikes and VTA 5-15 Hz oscillations. 

clear all; clear functions; clc

%%
addpath('/Users/parkj/Desktop/OldPC/Matlab');
addpath('/Users/parkj/Desktop/OldPC/Matlab/Functions');

%% load signal and timestamps modify this!
filedirectory = '/Users/parkj/Desktop/OldPC/Data/Anxiety_goal_directed_behavior/MAT/LFP/'; cd(filedirectory); % directory for evt timestamps, behavioral data (tstpbeh)
load('LFP_AGB_periispk_LFPtimeseries_PFCVTA_100715','tstpbeh'); % load behavioral time stamps        

filedirectory1 = '/Users/parkj/Desktop/OldPC/Data/Anxiety_goal_directed_behavior/MAT/'; cd(filedirectory1); % directory for evt timestamps, behavioral data (tstpbeh)
load('VTA_sua_dat_022115','dat'); % load LFP timeseries and behavioral time stamps        
VTAdat = dat; clearvars dat

filedirectory2 = '/Users/parkj/Desktop/OldPC/Data/Anxiety_goal_directed_behavior/NEX_LFP_Aligned/VTA'; cd(filedirectory2); % directory for .NEX files containing VTA LFP data
fileListVTAlfp = dir('*VTA*'); % VTA fileList

cd('/Users/parkj/Desktop/OldPC/Data/Anxiety_goal_directed_behavior/MAT/VTA_cell_classification') % directory where DA index is saved
load('VTA_SUA_unitClassification_120115','DAIdx')

cd(filedirectory2) % Back to the file directory

% get VTA file ID
for f = 1:length(fileListVTAlfp)
    VSpkLFPid{f,1} = strcat(fileListVTAlfp(f).name(1:10),fileListVTAlfp(f).name(16:17),fileListVTAlfp(f).name(36:39));  % the file ID being processed 
end
clearvars f 

%% set user variables for rasterdata modify these!
uv.binsize = 0.001;       % binsize in sec (1 msec, point process)
uv.win = [2.25 2.25];     % peri-event time window (-1 to 0 s with each event occurring at time = 0, 100 ms extra-period added to both ends for binning)
uv.edges = [-2.25:uv.binsize:2.249];   % 1 ms bins to be used for histc
uv.timelag = [-100:4:100]./1000;  % time lags in msec
uv.passband = [5 15];  % passband for filtering 
uv.detrend_win = [0.5, 0.05];    %detrending window size and step size 
uv.Fs = 1000;          % sampling frequency 

%% Time-lag Phase-locking analysis
VTAlfpidx = cell(length(VTAdat),1);     % cell containing the VTA LFP indices
PLVmat = nan(length(VTAdat),3,length(uv.timelag));    % Matrix containing the PLVs (unit x blocks (base, b1, b2, b3) x # of timelags)
PLVRaymat = nan(length(VTAdat),3,length(uv.timelag)); % Matrix containing the p-values from the Rayleigh's test for non-uniformity unit x blocks (base, b1, b2, b3) x # of timelags
basePLVmat = nan(length(VTAdat),3,length(uv.timelag));    % Matrix containing the baseline PLVs (unit x blocks (base, b1, b2, b3) x # of timelags)
basePLVRaymat = nan(length(VTAdat),3,length(uv.timelag)); % Matrix containing the baseline p-values from the Rayleigh's test for non-uniformity unit x blocks (base, b1, b2, b3) x # of timelags

for u = 1:length(VTAdat)    % the # of VTA units
    %% Get the peri-evt rasters 
    % Get the matching behavioral timestamps (tstpbeh)
    for o = 1:length(tstpbeh)
        behid(o,1) = strcmp(tstpbeh(o).name, VTAdat{u,1}(1:9));    % get the behavioral data index (logical)
    end
    clearvars o
    
    tempispkevt = tstpbeh(behid).ispk;    % get the current ispk timestamps using the tstpbehID
    tempcueevt = tstpbeh(behid).cue;      % get the current cue timestamps using the tstpbehID to be used for the baseline period in evtrawtimeserieswocutoff
    
    st = VTAdat{u,4}; % entire spike time stamps for the current unit 
    
    raster_data=zeros(length(uv.edges),length(tstpbeh(behid).ispk)); % preallocate the raster_data (trial x timebins, e.g. 150 x 4500)
    base_raster_data=zeros(length(uv.edges),length(tstpbeh(behid).ispk)); % preallocate the raster_data (trial x timebins, e.g. 150 x 4500) 
    for t = 1:length(tstpbeh(behid).ispk) % the number of trials
        [~, raster_data(:,t)] = perievtSpikeRaster( st, tstpbeh(behid).ispk(t,1), uv.win, uv.edges );    % raster_data:  e.g. 150x1000 matrix indicating spike (1) or not (0)  
        [~, base_raster_data(:,t)] = perievtSpikeRaster( st, tstpbeh(behid).ispk(t,1)-1, [4.5 0],[-4.5:uv.binsize:-0.001]);    % raster_data:  e.g. 150x1000 matrix indicating spike (1) or not (0)  
    end
    clearvars t   
    
    %% Get the matching VTA lfp data
    unitinfo = strcat(VTAdat{u,1}(1:10),VTAdat{u,1}(16)); % the file ID being processed
    unitch = str2num(VTAdat{u,3}(1:2)); % the ch# for the current unit
    
    for l = 1:length(VSpkLFPid) % the # of VTA LFP files
        tempidx(l,1) = strcmp(unitinfo,VSpkLFPid{l}(1:11)); % index for the VTA LFP matching to the current unit
    end
    clearvars l
    
    tempidx = logical(tempidx); % convert the index to logicals (index for the VTA LFP matching to the current unit)
    
    if sum(tempidx)>0 % In case, there's matching LFP file, proceed to the analysis
        
        if sum(tempidx)==1      % if there's only one matching LFP file (unilateral VTA)
            VTAlfpidx{u,1} = tempidx; % use it as the index
        elseif sum(tempidx)==2  % if there are two matching LFP files (bilateral VTA)
            if unitch <= 8      % left hemisphere
                tempidx(max(find(tempidx==1)))=0; % deselect the LFP data from right hemisphere
            elseif unitch > 8   % right hemisphere
                tempidx(min(find(tempidx==1)))=0; % deselect the LFP data from left hemisphere
            end
            VTAlfpidx{u,1} = tempidx; % use it as the index
        end
        
        % get the timeseries for the entire session
        [~,~,dts,fn,d] = nex_cont(fileListVTAlfp(tempidx).name, fileListVTAlfp(tempidx).name(36:39)); % get the continuous data of each animal
        [ crtd ] = crttimeseries(d, dts, fn);
        
        %% get PLVs across timelags  
        for g = 1:length(uv.timelag) % # of timelags (e.g. 51)
            % loop to get peri-evt (ISPK) psd using the function evtrawtimeseries containing the chronux function mtspectramc
            for t = 1:size(tempispkevt,1) % number of total trials 
                [ tempispklfp(:,t), tempbaselfp(:,t) ] = evtrawtimeserieswocutoff( crtd, tempispkevt(t,1) + uv.timelag(g), tempcueevt, -1 ); % -1 is the baseendpoint, baseline period will be the 4500 ms period during ITI ending 1 s before each ispk onset            
            end
            [PLVmat(u,:,g),PLVRaymat(u,:,g),pfPhasemat(u,:,g)] = newPLVrasterlfprandsample( tempispklfp, raster_data, uv, unitinfo(end)); % get the Phase-locking values
            [basePLVmat(u,:,g),basePLVRaymat(u,:,g),basepfPhasemat(u,:,g)] = newPLVrasterlfprandsample( tempbaselfp, base_raster_data, uv, unitinfo(end)); 
            fprintf('processed timelag #%d\n', g);     % report the progress (completion of each unit)  
        end
        clearvars t d dts fn
     
    else % In case, there's no matching LFP file, do not perform further analysis 
        
    end
            
    clearvars raster_data unit* temp*        
    fprintf('processed unit #%d\n', u);     % report the progress (completion of each unit)  
end
clearvars u

%% Post-processing
valcellidx = ~isnan(squeeze(sum(PLVmat,2))); % idx for val (non-NaN) cells
valcellidx = valcellidx(:,1); % idx for val (non-NaN) cells

valPLVmat = PLVmat(valcellidx,:,:); % periEVT PLVmat exclude NaNs
valbasePLVmat = basePLVmat(valcellidx,:,:); % baseline PLVmat exclude NaNs

% mean, SEM PLVs across time lags
[B1mPLV,~,B1sPLV] = meanstdsem(squeeze(valPLVmat(:,1,:)));
[B2mPLV,~,B2sPLV] = meanstdsem(squeeze(valPLVmat(:,2,:)));
[B3mPLV,~,B3sPLV] = meanstdsem(squeeze(valPLVmat(:,3,:)));

% get the proportion of significantly phase-locked neurons
valPLVRaymat = PLVRaymat(valcellidx,:,:);
sigproPLV = nan(size(valPLVRaymat,3),3); % matrix to contain the proportion of significantly phase-locked units 
sigPLV = nan(size(valPLVRaymat,3),3);    % matrix to contain the phase-locking values of the significant units
for lag = 1:size(valPLVRaymat,3) % # of time lags
    sigproPLV(lag,1) = sum(valPLVRaymat(:,1,lag) <= 0.05)/size(valPLVRaymat,1);
    sigproPLV(lag,2) = sum(valPLVRaymat(:,2,lag) <= 0.05)/size(valPLVRaymat,1);
    sigproPLV(lag,3) = sum(valPLVRaymat(:,3,lag) <= 0.05)/size(valPLVRaymat,1);
end
clearvars l

% significant phase-locking index
B1sigidx = valPLVRaymat(:,1,25)<0.05; % significant phase-locking at time 0
B2sigidx = valPLVRaymat(:,2,25)<0.05; % significant phase-locking at time 0
B3sigidx = valPLVRaymat(:,3,25)<0.05; % significant phase-locking at time 0

[B1sigmPLV,~,B1sigsPLV] = meanstdsem(squeeze(valPLVmat(B1sigidx,1,:))); % stats only including significant phase-locked units across time 
[B2sigmPLV,~,B2sigsPLV] = meanstdsem(squeeze(valPLVmat(B2sigidx,2,:))); % stats only including significant phase-locked units across time 
[B3sigmPLV,~,B3sigsPLV] = meanstdsem(squeeze(valPLVmat(B3sigidx,3,:))); % stats only including significant phase-locked units across time

[mB1fldincS,~,sB1fldincS] = meanstdsem(max(squeeze(valPLVmat(B1sigidx,1,:))./squeeze(valbasePLVmat(B1sigidx,1,:)),[],2)); % Fold increase from the baseline phase-locking level significantly phase-locked units
[mB2fldincS,~,sB2fldincS] = meanstdsem(max(squeeze(valPLVmat(B2sigidx,1,:))./squeeze(valbasePLVmat(B2sigidx,1,:)),[],2)); % Fold increase from the baseline phase-locking level significantly phase-locked units
[mB3fldincS,~,sB3fldincS] = meanstdsem(max(squeeze(valPLVmat(B3sigidx,1,:))./squeeze(valbasePLVmat(B3sigidx,1,:)),[],2)); % Fold increase from the baseline phase-locking level significantly phase-locked units

%B1sigmaxPLV = max(squeeze(valPLVmat(B1sigidx,1,:)),[],2);  % Block1 max PLVs of the significantly phase-locked units, peri-action period
%B1sigmaxPLVbase = max(squeeze(valbasePLVmat(B1sigidx,1,:)),[],2); % Block1 max PLVs of the significantly phase-locked units, baseline period

%B2sigmaxPLV = max(squeeze(valPLVmat(B2sigidx,1,:)),[],2);  % Block2 max PLVs of the significantly phase-locked units, peri-action period
%B2sigmaxPLVbase = max(squeeze(valbasePLVmat(B2sigidx,1,:)),[],2); % Block2 max PLVs of the significantly phase-locked units, baseline period

%B3sigmaxPLV = max(squeeze(valPLVmat(B3sigidx,1,:)),[],2);  % Block3 max PLVs of the significantly phase-locked units, peri-action period
%B3sigmaxPLVbase = max(squeeze(valbasePLVmat(B3sigidx,1,:)),[],2); % Block3 max PLVs of the significantly phase-locked units, baseline period

[mB1fldincNS,~,sB1fldincNS] = meanstdsem(max(squeeze(valPLVmat(~B1sigidx,1,:))./squeeze(valbasePLVmat(~B1sigidx,1,:)),[],2)); % Fold increase from the baseline phase-locking level
[mB2fldincNS,~,sB2fldincNS] = meanstdsem(max(squeeze(valPLVmat(~B2sigidx,1,:))./squeeze(valbasePLVmat(~B2sigidx,1,:)),[],2)); % Fold increase from the baseline phase-locking level
[mB3fldincNS,~,sB3fldincNS] = meanstdsem(max(squeeze(valPLVmat(~B3sigidx,1,:))./squeeze(valbasePLVmat(~B3sigidx,1,:)),[],2)); % Fold increase from the baseline phase-locking level

%B1NSmaxPLV = max(squeeze(valPLVmat(~B1sigidx,1,:)),[],2);  % Block1 max PLVs of the significantly phase-locked units, peri-action period
%B1NSmaxPLVbase = max(squeeze(valbasePLVmat(~B1sigidx,1,:)),[],2); % Block1 max PLVs of the significantly phase-locked units, baseline period

%B2NSmaxPLV = max(squeeze(valPLVmat(~B2sigidx,1,:)),[],2);  % Block2 max PLVs of the significantly phase-locked units, peri-action period
%B2NSmaxPLVbase = max(squeeze(valbasePLVmat(~B2sigidx,1,:)),[],2); % Block2 max PLVs of the significantly phase-locked units, baseline period

%B3NSmaxPLV = max(squeeze(valPLVmat(~B3sigidx,1,:)),[],2);  % Block3 max PLVs of the significantly phase-locked units, peri-action period
%B3NSmaxPLVbase = max(squeeze(valbasePLVmat(~B3sigidx,1,:)),[],2); % Block3 max PLVs of the significantly phase-locked units, baseline period

%% save
cd('/Users/parkj/Desktop/OldPC/Data/Anxiety_goal_directed_behavior/MAT/LFP/PLV')
saveName = 'PLVtimeLag_analysis_AGB_periISPK_VTAspk_VTAlfp_100ms_103015';
%save(saveName);
 
%% plot the time-lagged Phase-locking value
cmap(1,:)=[30 144 255]./255;  % B1 
cmap(2,:)=[0 238 118]./255;   % B2
cmap(3,:)=[220 20 60]./255;   % B3

% PLVs across time lags
figure;
boundedline(uv.timelag, smooth(B1mPLV,5), smooth(B1sPLV,5), uv.timelag, smooth(B2mPLV,5), smooth(B2sPLV,5), uv.timelag, smooth(B3mPLV,5), smooth(B3sPLV,5), 'alpha','cmap', cmap) 
ylim([0.09 0.11])
set(gca,'XTick',[-0.1:0.05:0.1],'FontSize',14);
set(gca,'TickDir','out')

%  figure;
%  boundedline(uv.timelag, smooth(B1sigmPLV,5), smooth(B1sigsPLV,5), uv.timelag, smooth(B2sigmPLV,5), smooth(B2sigsPLV,5), uv.timelag, smooth(B3sigmPLV,5), smooth(B3sigsPLV,5), 'alpha','cmap', cmap) 
%  ylim([0.1 0.15])
%  set(gca,'XTick',[-0.1:0.05:0.1],'FontSize',14);
%  set(gca,'TickDir','out')

% ur = resample(smooth(sigproPLV(:,1),5),1,2);

% Proportion Bars
figure;
hold on;
bar(uv.timelag, smooth(sigproPLV(:,1),5), 'FaceColor', cmap(1,:)); bar(uv.timelag, smooth(sigproPLV(:,3),5), 'FaceColor', cmap(3,:));    
xlim([-0.105 0.105])
set(gca,'XTick',[-0.1:0.05:0.1],'FontSize',14);
set(gca,'TickDir','out')

% Mean PLV bars
[B1premPLV,~,B1presPLV] = meanstdsem(mean(squeeze(valPLVmat(:,1,1:25)),2));
[B2premPLV,~,B2presPLV] = meanstdsem(mean(squeeze(valPLVmat(:,2,1:25)),2));
[B3premPLV,~,B3presPLV] = meanstdsem(mean(squeeze(valPLVmat(:,3,1:25)),2));
hold off;

figure;
hold on;
bar(1:3,[B1premPLV; B2premPLV; B3premPLV])
errorbar(1:3,[B1premPLV; B2premPLV; B3premPLV],[B1presPLV; B2presPLV; B3presPLV],'.')
ylim([0.09 0.11])
set(gca, 'YTick',[0.09:0.01:0.12])
ylabel('Phase-locking value')
hold off;
set(gca,'TickDir','out')

% % Mean PLV bars significantly phase-locked units
% [B1premPLVsig,~,B1presPLVsig] = meanstdsem(mean(squeeze(valPLVmat(B1sigidx,1,1:25)),2));
% [B2premPLVsig,~,B2presPLVsig] = meanstdsem(mean(squeeze(valPLVmat(B2sigidx,2,1:25)),2));
% [B3premPLVsig,~,B3presPLVsig] = meanstdsem(mean(squeeze(valPLVmat(B3sigidx,3,1:25)),2));
% hold off;
% 
% figure;
% hold on;
% bar(1:3,[B1premPLVsig; B2premPLVsig; B3premPLVsig])
% errorbar(1:3,[B1premPLVsig; B2premPLVsig; B3premPLVsig],[B1presPLVsig; B2presPLVsig; B3presPLVsig],'.')
% ylim([0.09 0.14])
% set(gca, 'YTick',[0.09:0.01:0.15])
% ylabel('Phase-locking value')
% set(gca, 'TickDir','out')
% hold off;

% Mean PLV bars Phase-locked neurons only
[B1premPLVsig,~,B1presPLVsig] = meanstdsem(mean(squeeze(valPLVmat(B1sigidx,1,1:25)),2));
[B2premPLVsig,~,B2presPLVsig] = meanstdsem(mean(squeeze(valPLVmat(B2sigidx,2,1:25)),2));
[B3premPLVsig,~,B3presPLVsig] = meanstdsem(mean(squeeze(valPLVmat(B3sigidx,3,1:25)),2));
hold off;

hold on;
bar(1:3,[B1premPLVsig; B2premPLVsig; B3premPLVsig])
errorbar(1:3,[B1premPLVsig; B2premPLVsig; B3premPLVsig],[B1presPLVsig; B2presPLVsig; B3presPLVsig],'.')
ylim([0.09 0.14])
set(gca, 'YTick',[0.09:0.01:0.15])
ylabel('Phase-locking value')
hold off;
set(gca,'TickDir','out')

% Mean PLV fold increase from the baseline level
figure;
hold on;
bar(1:6,[mB1fldincNS; mB2fldincNS; mB3fldincNS; mB1fldincS; mB2fldincS; mB3fldincS]) 
errorbar(1:6,[mB1fldincNS; mB2fldincNS; mB3fldincNS; mB1fldincS; mB2fldincS; mB3fldincS],[sB1fldincNS; sB2fldincNS; sB3fldincNS; sB1fldincS; sB2fldincS; sB3fldincS],'.') 
ylim([0.5 2])
hold off;
set(gca, 'YTick',[0.5:0.5:2])
set(gca,'TickDir','out')

% Colormap indicating the scaled max and min PLVs across time lags (including only valid significantly phase-locked units)
B1valsigPLVlag = squeeze(valPLVmat(B1sigidx,1,:)); % B1 PLV val significant units across time lags
B2valsigPLVlag = squeeze(valPLVmat(B2sigidx & ~B1sigidx,2,:)); % B2 PLV val significant units across time lags
B3valsigPLVlag = squeeze(valPLVmat(B3sigidx & ~(B1sigidx | B2sigidx),3,:)); % B3 PLV val significant units across time lags

valsigPLVlag = [B1valsigPLVlag; B2valsigPLVlag; B3valsigPLVlag];

maxvalsigPLVlag = max(valsigPLVlag,[],2); % max B1valsigPLVlag across time lags
minvalsigPLVlag = min(valsigPLVlag,[],2); % min B1valsigPLVlag across time lags
maxminPLV = (valsigPLVlag-repmat(minvalsigPLVlag,[1 size(valsigPLVlag,2)]))./repmat((maxvalsigPLVlag-minvalsigPLVlag),[1 size(valsigPLVlag,2)]);
[~,maxidx] = max(maxminPLV,[],2); % find the time lag with the max PLV value
% maxidx-27 % time bins with the max PLV value
srtmaxminPLV = sortrows([maxminPLV,maxidx],size(valsigPLVlag,2)+1); % sortrows based on the location of the max time lag in an ascending manner
for u = 1:size(srtmaxminPLV,1) % increment valid significantly phase-locked units
    smsrtmaxminPLV(u,:) = smooth(srtmaxminPLV(u,1:end-1),5); % smooth
end

imagesc(smsrtmaxminPLV(:,1:end-1))
colormap jet
set(gca,'TickDir','out')

%% Proportion of DA and non-DA cells phase-locked to the local 10 Hz oscillation
cd('/Users/parkj/Desktop/OldPC/Data/Anxiety_goal_directed_behavior/MAT/VTA_cell_classification') % directory where DA index is saved
load('VTA_SUA_unitClassification_120115','DAIdx')
cd(filedirectory2) % Back to the file directory

valDAIdx = DAIdx(valcellidx); % valid DA & non-DA index
valDApl = B1sigidx & valDAIdx; % valid DA cells phase-locked
valnDApl = B1sigidx & ~valDAIdx; % valid DA cells phase-locked

valDAplper = sum(valDApl)/sum(valDAIdx); % percent of phase-locked DA units
valnDAplper = sum(valnDApl)/sum(~valDAIdx); % percent of phase-locked non-DA units

%% Individual animal VTA plv comparison
VTAmID = zeros(length(VTAdat),1);   % VTA animal ID 
VTAmID(1:28,1)   = 0;  % AD01_Day1
VTAmID(29:48,1)  = 1;  % AD01_Day3
VTAmID(49:56,1)  = 2;  % AD02_Day4
VTAmID(57:60,1)  = 3;  % AD03_Day3
VTAmID(61:67,1)  = 5;  % AD05_Day2
VTAmID(68:71,1)  = 6;  % AD06_Day2
VTAmID(72:79,1)  = 11; % JP11_Day2
VTAmID(80:89,1)  = 12; % JP12_Day4
VTAmID(90:102,1) = 14; % JP14_Day5

VTAmbymPLV      = cell(length(unique(VTAmID)),1); % VTA animal-by-animal PLV 
VTAmbymPLVsig   = cell(length(unique(VTAmID)),1); % VTA animal-by-animal PLV significant animals only 
uniqM           = unique(VTAmID);   % the list of unique animal ID

countM        = 1;  % animal counter - to count the number of animals
countUnitPerM = 1;  % unit counter per animal - to count the number of units per animal

for u = 1:length(PLVmat) 
    
    if u == 1 % in case of the 1st unit
       countM = 1;        % set the animal number to 1
       countUnitPerM = 1; % set the unit number to 1
       
       VTAmbymPLV{countM}(countUnitPerM,1) = u;        % put unit number
       VTAmbymPLV{countM}(countUnitPerM,2) = DAIdx(u); % put DA index
       VTAmbymPLV{countM}(countUnitPerM,3) = PLVRaymatIdx(u); % put PLVRaymatIdx index 
       VTAmbymPLV{countM}(countUnitPerM,4) = max(PLVmat(u,1,:)); % block 1 PLV max value across different lags
       VTAmbymPLV{countM}(countUnitPerM,5) = max(PLVmat(u,2,:)); % block 2 PLV max value across different lags
       VTAmbymPLV{countM}(countUnitPerM,6) = max(PLVmat(u,3,:)); % block 3 PLV max value across different lags
       
       
    elseif VTAmID(u,1)==VTAmID(u-1,1) % in case the animal number stays the same
       %countM = countM;  % animal count does NOT increment  
       countUnitPerM = countUnitPerM + 1; % unit count per M increment
       
       VTAmbymPLV{countM}(countUnitPerM,1) = u;        % put unit number
       VTAmbymPLV{countM}(countUnitPerM,2) = DAIdx(u); % put DA index
       VTAmbymPLV{countM}(countUnitPerM,3) = PLVRaymatIdx(u); % put PLVRaymatIdx index
       VTAmbymPLV{countM}(countUnitPerM,4) = max(PLVmat(u,1,:)); % block 1 PLV max value across different lags
       VTAmbymPLV{countM}(countUnitPerM,5) = max(PLVmat(u,2,:)); % block 2 PLV max value across different lags
       VTAmbymPLV{countM}(countUnitPerM,6) = max(PLVmat(u,3,:)); % block 3 PLV max value across different lags
        
    elseif VTAmID(u,1)~=VTAmID(u-1,1) % in case the animal number of the current unit differs from the previous unit 
        countM = countM + 1; % must increment the M # 
        countUnitPerM = 1;   % reset the # of unit to 1
       
       VTAmbymPLV{countM}(countUnitPerM,1) = u;        % put unit number
       VTAmbymPLV{countM}(countUnitPerM,2) = DAIdx(u); % put DA index
       VTAmbymPLV{countM}(countUnitPerM,3) = PLVRaymatIdx(u); % put PLVRaymatIdx index
       VTAmbymPLV{countM}(countUnitPerM,4) = max(PLVmat(u,1,:)); % block 1 PLV max value across different lags
       VTAmbymPLV{countM}(countUnitPerM,5) = max(PLVmat(u,2,:)); % block 2 PLV max value across different lags
       VTAmbymPLV{countM}(countUnitPerM,6) = max(PLVmat(u,3,:)); % block 3 PLV max value across different lags 
        
    end
         
end

allMplvDA = nan(length(VTAmbymPLV),1); % PLVs of all animals
sigMplvDA = nan(length(VTAmbymPLV),1); % PLVs of animals with significant PLVs

for m = 1:length(VTAmbymPLV) % increment animals
  
    if m ~= 5
        allMplvDA(m,1) = nanmean(max(VTAmbymPLV{m,1}(logical(VTAmbymPLV{m,1}(:,2)),4:6),[],2)); % mean plv of all DA neurons animal-by-animal
        sigMplvDA(m,1) = nanmean(max(VTAmbymPLV{m,1}(logical(VTAmbymPLV{m,1}(:,2))&logical(VTAmbymPLV{m,1}(:,3)),4:6),[],2)); % mean plv of signifiantly phase-locked DA neurons animal-by-animal
    elseif m == 5
        allMplvDA(m,1) = nanmean(nanmedian(VTAmbymPLV{m,1}(logical(VTAmbymPLV{m,1}(:,2)),4:6),2)); % mean plv of all DA neurons animal-by-animal
        sigMplvDA(m,1) = nanmean(nanmedian(VTAmbymPLV{m,1}(logical(VTAmbymPLV{m,1}(:,2))&logical(VTAmbymPLV{m,1}(:,3)),4:6),2)); % mean plv of signifiantly phase-locked DA neurons animal-by-animal
    end
    
    hold on;
    
    if m ~= 5
        for u = 1:size(VTAmbymPLV{m,1},1)
            if logical(VTAmbymPLV{m,1}(u,2))
                plot(m,max(VTAmbymPLV{m,1}(u,4:6),[],2),'-o','Color',[117 172 66]./255,'MarkerFaceColor',[117 172 66]./255,'MarkerSize',8)
            end
        end
    elseif m == 5
        for u = 1:size(VTAmbymPLV{m,1},1)
            if logical(VTAmbymPLV{m,1}(u,2))
                plot(m,max(VTAmbymPLV{m,1}(u,4:6),[],2),'-o','Color',[117 172 66]./255,'MarkerFaceColor',[117 172 66]./255,'MarkerSize',8)
            end
        end
    end
end

set(gca,'tickDir','out');
xlim([0 10])
ylim([0.08 0.3])
set(gca,'ytick',0:0.05:0.3);   
















