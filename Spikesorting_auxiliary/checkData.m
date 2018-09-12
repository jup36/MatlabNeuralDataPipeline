% I want to compare the source clusters and clustered units.
% How should I do it... Ahhhhhhh

clc; clearvars; close all;


%% Variable
source_file = 'E:\simData.mat';
jrc_file = 'E:\simData_jrc_20180911\simData.ap_imec3_opt3_jrc.mat';
phy_dir = 'E:\simData_kilosort_20180911\';


% imec
channel = 1:384;
ref_sites = [37 76 113 152 189 228 265 304 341 380];
channel(ref_sites) = []; 


%% loading source
load(source_file);
nC = length(spikeTime);
nSpk = cellfun(@length, spikeTime);

timeSrc = cell2mat(spikeTime);
cluSrc = cell2mat(cellfun(@(x, y) ones(length(x), 1)*y, spikeTime, num2cell(cellChannel'), 'UniformOutput', false));

[timeSrc, sortIdx] = sort(timeSrc);
cluSrc = cluSrc(sortIdx);
wavSrc = squeeze(wavSample(:, 1, :));
% source output: spikeTime, timeSrc, cluSrc, cellChannel, peakAmplitude,
% wavSample, nC

%% loading jrc
load(jrc_file, 'S_clu', 'spikeTimes');
inJrc = ~strcmp(S_clu.clusterNotes, 'noise'); % I excluded unsatisfactory units by tagging 'noise'.

spikeJrcTemp = cell(S_clu.nClusters, 1);
for iC = 1:S_clu.nClusters
    spikeJrcTemp{iC} = spikeTimes(S_clu.spikeClusters==iC);
end
spikeTimeJrc = spikeJrcTemp(inJrc);
cellChannelJrc = channel(S_clu.clusterSites(inJrc));
wavJrc = S_clu.trWav_raw_clu(:, :, inJrc);

timeJrc = cell2mat(spikeTimeJrc);
cluJrc = cell2mat(cellfun(@(x, y) ones(length(x), 1)*y, spikeTimeJrc, num2cell(cellChannelJrc'), 'UniformOutput', false));

[timeJrc, sortIdx] = sort(timeJrc);
cluJrc = cluJrc(sortIdx);
% jrc output: spikeTimeJrc, timeJrc, cluJrc, cellChannelJrc, wavJrc, nJ


%% kilosort-phy (you need https://github.com/kwikteam/npy-matlab)
spike_times = readNPY(fullfile(phy_dir, 'spike_times.npy')); % timestamp of all spikes, (n_spike, 1)
spike_template = readNPY(fullfile(phy_dir, 'spike_templates.npy')); % automated cluster of all spikes, (n_spike, 1)
spike_clusters = readNPY(fullfile(phy_dir, 'spike_clusters.npy')); % final cluster of all spikes, (n_spike, 1)
template = readNPY(fullfile(phy_dir, 'templates.npy')); % template waveform of all clusters, (n_original_cluster, n_timepoint, all_valid_channel)

% I find the main channel of each template which has the lowest amplitude.
[~, mainSiteTemplate] = min(min(template, [], 2), [], 3);

% load cluster data...
cn_name = fullfile(phy_dir, 'cluster_group.tsv');
fid = fopen(cn_name, 'r');
cn = textscan(fid, '%f%s%[^\n\r]', 'Delimiter', '\t', 'TextType', 'string', 'Headerlines', 1, 'EndOfLine', '\r\n');
fclose(fid);

% pick only good units
inPhy = strcmp(cn{2}, 'good');
unitNumber = cn{1}(inPhy); % final good unit number list
nP = length(unitNumber);

cellChannelPhy = zeros(nP, 1);
wavPhy = zeros(nP, 61, 'single');
timePhy = cell(nP, 1);
cluPhy = cell(nP, 1);
for iP = 1:nP
    spikeIndex = spike_clusters == unitNumber(iP);
    
    % final cluster can be composed of multiple templates. I picked the
    % most frequent template.
    spikeMainTemplateIndex = mode(spike_template(spikeIndex)) + 1; % spike_template start from 0
    
    % convert to imec channel
    cellChannelPhy(iP) = channel(mainSiteTemplate(spikeMainTemplateIndex));
    
    % get only the main channel waveform...
    wavPhy(iP, :) = squeeze(template(spikeMainTemplateIndex, 22:82, mainSiteTemplate(spikeMainTemplateIndex)));
    
    timePhy{iP} = double(spike_times(spikeIndex));
    cluPhy{iP} = ones(sum(spikeIndex), 1) * cellChannelPhy(iP);
end
spikeTimePhy = timePhy;
timePhy = cell2mat(timePhy);
cluPhy = cell2mat(cluPhy);

[timePhy, sortIdx] = sort(timePhy);
cluPhy = cluPhy(sortIdx);
% phy output: timePhy, cluPhy, cellChannelPhy, wavPhy, nP

clearvars -except spikeTime timeSrc cluSrc cellChannel peakAmplitude wavSrc nC nSpk...
    spikeTimeJrc timeJrc cluJrc cellChannelJrc wavJrc nJ ...
    spikeTimePhy timePhy cluPhy cellChannelPhy wavPhy nP
%% Analyze Phy
% source output: timeSrc, cluSrc, cellChannel, peakAmplitude, wavSample
% jrc output: timeJrc, cluJrc, cellChannelJrc, wavJrc
% phy output: timePhy, cluPhy, cellChannelPhy, wavPhy

matchPhy = cell(nC, 10);
rhoPhy = cell(nC, 1);
perfPhy = NaN(nC, 10);
for iC = 1:nC
    inSite = find(abs(cellChannelPhy - cellChannel(iC))<=4);
    nIn = length(inSite);
    
    rhoPhy{iC} = zeros(nIn, 2);
    rhoPhy{iC}(:, 1) = inSite;
    
    iIndex = 1;
    for iI = 1:nIn
        rhoTemp = zeros(1, 10);
        for iR = 1:10
            rhoTemp = corr(wavSrc(10:61, iC), wavPhy(inSite(iI), (iR:(iR+51)))');
        end
        rhoPhy{iC}(iI, 2) = max(rhoTemp);
        
        if rhoPhy{iC}(iI) >= 0.7
            matchPhy{iC, iIndex} = zeros(nSpk(iC), 4, 'single');
            
            for iS = 1:nSpk(iC)
                [minValue, minIndex] = min(abs(spikeTimePhy{inSite(iI)} - spikeTime{iC}(iS)));
                if minValue <= 15
                    matchPhy{iC, iIndex}(iS, :) = [inSite(iI), minIndex, spikeTime{iC}(iS), spikeTimePhy{inSite(iI)}(minIndex)];
                end
            end
            perfPhy(iC, iIndex) = mean(matchPhy{iC, iIndex}(:, 1)>0);
            iIndex = iIndex + 1;
        end
    end
end

[performancePhy, phyIndex] = max(perfPhy, [], 2);


%% Analyze JRC
matchJrc = cell(nC, 10);
rhoJrc = cell(nC, 1);
perfJrc = NaN(nC, 10);
for iC = 1:nC
    inSite = find(abs(cellChannelJrc - cellChannel(iC))<=4);
    nIn = length(inSite);
    
    rhoJrc{iC} = zeros(nIn, 2);
    rhoJrc{iC}(:, 1) = inSite;
    
    iIndex = 1;
    for iI = 1:nIn
        rhoTemp = zeros(1, 10);
        for iR = 1:10
            rhoTemp = corr(wavSrc(10:61, iC), wavJrc(iR:(iR+51), 1, inSite(iI)));
        end
        rhoJrc{iC}(iI, 2) = max(rhoTemp);
        
        if rhoJrc{iC}(iI) >= 0.7
            matchJrc{iC, iIndex} = zeros(nSpk(iC), 4, 'single');
            
            for iS = 1:nSpk(iC)
                [minValue, minIndex] = min(abs(spikeTimeJrc{inSite(iI)} - spikeTime{iC}(iS)));
                if minValue <= 15
                    matchJrc{iC, iIndex}(iS, :) = [inSite(iI), minIndex, spikeTime{iC}(iS), spikeTimeJrc{inSite(iI)}(minIndex)];
                end
            end
            perfJrc(iC, iIndex) = mean(matchJrc{iC, iIndex}(:, 1)>0);
            iIndex = iIndex + 1;
        end
    end
end

[performanceJrc, jrcIndex] = max(perfJrc, [], 2);


%% check false positive
nPhy = cellfun(@length, spikeTimePhy);
nJrc = cellfun(@length, spikeTimeJrc);
falsePhy = cell(nP, 1);
falseJrc = cell(nJ, 1);
for iC = 1:nC
    inCluPhy = find(matchPhy{iC, phyIndex(iC)}(:, 1)>0);
    if ~isempty(inCluPhy)
        targetClusterPhy = matchPhy{iC, phyIndex(iC)}(inCluPhy(1), 1);
        nFound = sum(matchPhy{iC, phyIndex(iC)}(:, 1)>0);
        falsePhy{targetClusterPhy} = [falsePhy{targetClusterPhy}, 1 - (nFound / nPhy(targetClusterPhy))];
    end
    
    inCluJrc = find(matchJrc{iC, jrcIndex(iC)}(:, 1)>0);
    if ~isempty(inCluJrc)
        targetClusterJrc = matchJrc{iC, jrcIndex(iC)}(inCluJrc(1), 1);
        nFound = sum(matchJrc{iC, jrcIndex(iC)}(:, 1)>0);
        falseJrc{targetClusterJrc} = [falseJrc{targetClusterJrc}, 1 - (nFound / nJrc(targetClusterJrc))];
    end
end

falsePhy(cellfun(@isempty, falsePhy)) = {NaN};
falseJrc(cellfun(@isempty, falseJrc)) = {NaN};


subplot(2, 2, 1);
plot(-peakAmplitude, performancePhy, '.');
title('Kilosort performance');

subplot(2, 2, 2);
plot(-peakAmplitude, performanceJrc, '.');
title('Jrc performance');

subplot(2, 2, 3);
plot(-min(wavPhy, [], 2), cellfun(@min, falsePhy), '.'); % I don't know what's the phy's waveform unit is but it is very small.
title('Kilosort false positive / overmerged');

subplot(2, 2, 4);
plot(-min(squeeze(wavJrc(:, 1, :)), [], 1)', cellfun(@min, falseJrc), '.');
title('Jrc false positive / overmerged');
