addpath(genpath('E:\MatlabNeuralDataPipeline'))
%filePath = '/Volumes/RAID2/parkj/NeuralData/ITphys/IT01_Ldms_M1_121317/Matfiles';
cd(filePath)

binFile = dir(fullfile(filePath,'*.lf.bin')); % look for nidq.bin file   
binName = binFile.name;

meta = ReadMeta(binName, filePath); % get the meta data (structure)

% load ts
load('BehVariables','ts')

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
uv.binsize = 0.001;     % binsize in sec (1 msec, point process)
uv.win = [2.25 2.25];   % peri-event time window (-1 to 0 s with each event occurring at time = 0, 100 ms extra-period added to both ends for binning)
uv.edges = [-2.25:uv.binsize:2.249];    % 1 ms bins to be used for histc
uv.passband = [5 15];  %passband for filtering 
uv.detrend_win = [0.5, 0.05];    %detrending window size and step size 
uv.Fs = 1000;          %sampling frequency 

% Read the binary data (entire samples)
nSamp         = SampRate(meta);          % sampling rate (default: 25kHz) 
totalTimeSecs = str2double(meta.fileTimeSecs); % total duration of file in seconds

lfpCh = [150, 330];
lfpData = zeros(2 ,floor(totalTimeSecs*1000)); %zeros(str2double(meta.nSavedChans)-1 ,floor(totalTimeSecs*1000)); % the time resolution will be 1000Hz (1ms) after decimation

for t = 0:totalTimeSecs-1 % read second-by-second incrementally to avoid a memory issue
    tempDataArray = ReadBin(t*nSamp, nSamp, meta, binName, filePath); % read bin data for each second
    
    for ch = 1:length(lfpCh) %size(lfpData,1) 
        resampTempDataArray = resample(tempDataArray(ch,:),10,1); % resample as 25kHz (the default sample rate for LFP data is 2.5kHz)
        lfpData(ch,t*1000+1:(t+1)*1000) = decimate(resampTempDataArray,round(nSamp*10/1000));
    
    end
    fprintf('processed %d\n', t+1) 
end
clearvars i

% Gain correction for channnels of interest
lfpSTR = GainCorrectIM(lfpData(1,:), 1, meta);   % gain-corrected lfp 
lfpCTX = GainCorrectIM(lfpData(2,:), 1, meta);   

%% get peri-event spectral densities and normalize them relative to the baseline psdS
for i = 1:size(lfpData,1)
    meand = mean(lfpData(i,:),2);      % get the mean to set the cut-off criteria
    stdd = std(lfpData(i,:),0,2);      % get the std to set the cut-off criteria
    upperlimit = meand + 5*stdd;
    lowerlimit = meand - 5*stdd;
    
    tempbase = cell(1,1);   % tempbase contains valid baselines  
    tempbasecount = 0;      % counter for valid tempbase
          
    %% loop to get peri-evt (reach) psdS using the function longeventpsdS containing the chronux function mtspectramc    
    for ii = 1:length(ts.reachStart) % number of total reachs - 1 (# of trials), because the last trial was most likely to be truncated 
        [ psdS.reach{ii,1,i}, raw.reach{ii,1,i}, psdS.norm_reach{ii,1,i}, psdS.base_reach{ii,1,i}, raw.base_reach{ii,1,i}, t, f ] = eventpsdAGB( lfpData(i,:), ts.reachStart(ii)./1000, ts.reward./1000+2, movingwin, params, sv, upperlimit, lowerlimit );       % get t and f to use it for later in the script
        if sum(isnan(psdS.base_reach{ii,1,i}))==0
           tempbasecount = tempbasecount + 1;
           tempbase{tempbasecount,1}(:,:) = psdS.base_reach{ii,1,i};
        else 
        end
     
    end
    clearvars ii j
    
    % normalization to the avg collected baseline
    [psdS.normavg_reach(1:length(ts.reachStart),1,i), collectbasemat(i,:), collectbasestdmat(i,:)] = psdnormalization(psdS.reach(1:length(ts.reachStart),:,i), tempbase, t, f);       % number of total reachs - 1 (# of trials), because the last trial was most likely to be truncated 
    
end
clearvars i

%% image the normalized psd
valTrialCount = 0;
for i = 1:size(lfpData,1)
    for ii = 1:size(psdS.normavg_reach,1)
        if ~isnan(psdS.normavg_reach{ii,1,i})
            valTrialCount = valTrialCount + 1;
            psdS.collectNormAvgReach{1,i}(:,:,valTrialCount) = psdS.normavg_reach{ii,1,i};
            psdS.collectRawReachLFP{1,i}(:,valTrialCount) = raw.reach{ii,1,i}';
        end
    end
    valTrialCount = 0;
end

%%
lfpdtSTR = locdetrend(psdS.collectRawReachLFP{1,1},uv.Fs, uv.detrend_win);
lfpdtCTX = locdetrend(psdS.collectRawReachLFP{1,2},uv.Fs, uv.detrend_win);
t = (1:size(lfpdtSTR, 1))'/uv.Fs; % time axis
filterlfpSTR = FilterLFP([t, lfpdtSTR], 'passband', uv.passband);
filterlfpCTX = FilterLFP([t, lfpdtCTX], 'passband', uv.passband);
        % Filter VTA lfp
        filter_lfp = FilterLFP([t, lfpdt], 'passband', uv.passband);
        % Instantaneous phase
        SpkLFP(f).PispklfpPhase = nan(length(uv.edges),size(SpkLFP(f).Pispklfp,2)); % matrix to contain the instantaneous phases
        SpkLFP(f).PispklfpPhase = angle(hilbert(filter_lfp(:, 2:end))); % instantaneous phase of oscillation

%% select files to be plotted
plotvalpsdmat = nanmean(psdS.collectNormAvgReach{2},3);    % average across rats 
%plotvalpsdmat = nanmean(eachvalpsdcell{25},3);

%% smoothing for the averaged psd in both time and frequency dimensions
for i = 1:size(plotvalpsdmat,2)         % for smoothing in time domain
    sm_avgvalpsdmat(:,i) = smooth(plotvalpsdmat(:,i),5);      % smoothing across row (time)
end
clearvars i

%% image using surf figure
[n_row, n_col]=size(plotvalpsdmat);
time_bin=4;      % 4 sec, the length of time window     
max_freq=100;

time_scale = (0+time_bin/n_row)-2:time_bin/n_row:(time_bin)-2;             % the length of the entire timewindow is 4 sec, which is divided into 2 + 2 around the event. 
 
[x,y]=meshgrid(time_scale,[max_freq/n_col:max_freq/n_col:max_freq]);
 
x = x';         % to match dimension
y = y';         % to match dimension

figure;
figure = surf(x,y,sm_avgvalpsdmat);
set(figure,'edgecolor','none')
view(0,90);

set(gca,'xtick',-1:1:1);
set(gca,'ytick',10:10:100);
set(gca,'TickDir','out');
set(gca,'FontSize',14);
%set(gca,'Layer','none');
set(gca,'GridLineStyle','none');
colorbar('EastOutside');
caxis([-1 1.5])
axis tight

xlabel('Time (sec)','FontSize',14)
ylabel('Frequency (Hz)','FontSize',14)   





