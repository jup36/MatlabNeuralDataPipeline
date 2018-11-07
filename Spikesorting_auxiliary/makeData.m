clc; clearvars; close all;

%% Variables
% recording
gain = 500;
nCell = 200;
spikeVariability =  0.05; % percent
recordingDuration = 10*60; % second
noiseAmplitude = 6; % microvolts
frRange = [0.1, 10];
samplingRate = 30000;


% voltage
bitPerMicrovolt = 2^10 / (1.2 / gain * 1000000); % 2.3438 uV/bit = 0.4267 bit/uV


% saving directory
saveDir = 'E:';

%% IMEC probe information
nChannel = 384;
channel = 1:nChannel;
geometry = zeros(nChannel, 2);
viHalf = 0:(nChannel/2-1);
geometry(1:2:end,2) = viHalf * 20;
geometry(2:2:end,2) = geometry(1:2:end,2);
geometry(1:4:end,1) = 16; %0 before
geometry(2:4:end,1) = 48; %32 before
geometry(3:4:end,1) = 0;  %16 before
geometry(4:4:end,1) = 32; %48 before
ref_sites = [37 76 113 152 189 228 265 304 341 380];
channel(ref_sites) = []; 
geometry(ref_sites,:) = [];

nLiveChannel = length(channel);
nearestChannel = zeros(nLiveChannel, 14);
for iC = 1:nLiveChannel
    distance = (geometry(:, 1) - geometry(iC, 1)).^2 + (geometry(:, 2) - geometry(iC, 2)).^2;
    [~, nearIndex] = sort(distance);
    nearestChannel(iC, :) = nearIndex(1:14);
end

% output: nChannel, nLiveChannel, nearestChannel

%% Spike time
fr = rand(nCell, 1) * (frRange(2) - frRange(1)) + frRange(1);
cellSite = randperm(nLiveChannel, nCell);
cellChannel = channel(cellSite);

spikeTime = cell(nCell, 1);
spikeIndex = cell(nCell, 1);
for jC = 1:nCell
    iti = round(geornd(fr(jC)/samplingRate, ceil(2*fr(jC)*recordingDuration), 1));
    iti(iti < 0.002*samplingRate) = [];
    spikeTemp = cumsum(iti);
    spikeTemp(spikeTemp > recordingDuration * samplingRate) = [];
    clipped = abs(spikeTemp/samplingRate - round(spikeTemp / samplingRate)) <= (64 / samplingRate);
    spikeTemp(clipped) = [];
    
    nSpike = length(spikeTemp);
    spikeTime{jC} = spikeTemp;
    spikeIndex{jC} = ones(nSpike, 1) * jC;
end

spikeAll = [cell2mat(spikeTime), cell2mat(spikeIndex)];
spikeAll = sortrows(spikeAll);

spikeBinned = cell(recordingDuration, 1);
for iT = 1:recordingDuration
    si = (spikeAll(:, 1) > (iT-1)*samplingRate) & (spikeAll(:, 1) < iT*samplingRate);
    spikeBinned{iT, 1} = spikeAll(si, :);
end
disp('spike time generated');
% output: spikeTime, nSpike, spikeBinned, cellChannel

%% Load spike waveform
load('wavSample.mat', 'wavSample');
cellIndex = randperm(size(wavSample, 3), nCell);
wavSample = wavSample(:, :, cellIndex);

% normalize waveform so that each waveform starts and ends at 0 uV.
% if it is not starting and ending from 0 uV, JRClust will detect that
% change and cluster it is a spike.
for iC = 1:nCell
    intercept = wavSample(1, :, iC);
    term = wavSample(end, :, iC);
    slopePerSample = (term - intercept) / 62;
    for iS = 1:63
        wavSample(iS, :, iC) = (wavSample(iS, :, iC) - intercept) - slopePerSample * (iS - 1);
    end
end

peakAmplitude = squeeze(min(wavSample(:, 1, :)));

disp('waveform loaded');
% output: wavSample

%save(fullfile(saveDir, 'simData.mat'), 'spikeTime', 'cellChannel', 'wavSample', 'peakAmplitude');

%% Generate data
disp('starting data generation');
fid = fopen(fullfile(saveDir, 'simData.bin'), 'w');
for iT = 1:recordingDuration
    data = randn(nChannel + 1, samplingRate, 'single') * (noiseAmplitude * bitPerMicrovolt);
    data(ref_sites, :) = 0;
    
    for iC = 1:length(spikeBinned{iT})
        activeChannel = channel(nearestChannel(cellSite(spikeBinned{iT}(iC, 2)), :));
        activeTime = spikeBinned{iT}(iC, 1) + (-15:47) - (iT-1)*samplingRate;
        spikeAmplitude = (1 + randn() * spikeVariability) * bitPerMicrovolt;
        
        data(activeChannel, activeTime) = data(activeChannel, activeTime) + spikeAmplitude * squeeze(wavSample(:, :, spikeBinned{iT}(iC, 2)))';
    end
    
    %dataWrite = int16(data(:));
    %fwrite(fid, dataWrite, 'int16');    
end
%fclose(fid);

% check data
hold on;
for iC = 1:100
    plot((1:3000)/30000, iC*100 + data(iC, 1:3000), 'Color', [0.5 0.5 0.5]);
end

figure();
for iS = 1:14
    subplot(14, 1, iS);
    hold on;
    plot(wavSample(:, iS, 1));
    plot([1, 63], [0, 0]);
end
