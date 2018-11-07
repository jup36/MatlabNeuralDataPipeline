% I want to compare the source clusters and clustered units.
% How should I do it... Ahhhhhhh

clc; clearvars; close all;

%% Variable
source_file = 'C:\Data\Spikes\SyntheticData_Dohoung\simData.mat';
jrc_file = 'C:\Data\Spikes\SyntheticData_Dohoung_jrcRun\simData_jrcRun.ap_imec3_opt3_jrc.mat'; % directory for jrc results
phy_dir = 'C:\Data\Spikes\SyntheticData_Dohoung'; % directory for kilosort results

% imec
channel = 1:384;
ref_sites = [37 76 113 152 189 228 265 304 341 380];
channel(ref_sites) = []; 

%% loading source
load(source_file);
nC = length(spikeTime); % # of cells
nSpk = cellfun(@length, spikeTime);

timeSrc = cell2mat(spikeTime);
cluSrc = cell2mat(cellfun(@(x, y) ones(length(x), 1)*y, spikeTime, num2cell(cellChannel'), 'UniformOutput', false));

[timeSrc, sortIdx] = sort(timeSrc);
cluSrc = cluSrc(sortIdx);
wavSrc = squeeze(wavSample(:, 1, :));
frSrc = nSpk./600; 
% source output: spikeTime, timeSrc, cluSrc, cellChannel, peakAmplitude,
% wavSample, nC

%% loading jrc
load(jrc_file, 'S_clu', 'viTime_spk');
jrcSpkTimes = viTime_spk; clearvars viTime_spk
%inJrc = ~strcmp(S_clu.clusterNotes, 'noise'); % I excluded unsatisfactory units by tagging 'noise'.

spikeJrcTemp = cell(S_clu.nClu, 1);
for iC = 1:S_clu.nClu
    spikeJrcTemp{iC} = jrcSpkTimes(S_clu.viClu==iC);
end
spikeTimeJrc = spikeJrcTemp(:);
cellChannelJrc = channel(S_clu.viSite_clu(:));
wavJrc = S_clu.trWav_raw_clu(:, :, :);

timeJrc = cell2mat(spikeTimeJrc);
cluJrc = cell2mat(cellfun(@(x, y) ones(length(x), 1)*y, spikeTimeJrc, num2cell(cellChannelJrc'), 'UniformOutput', false));

[timeJrc, sortIdx] = sort(timeJrc);
cluJrc = cluJrc(sortIdx); % the channel from 1:384 corresponding to the cluster 
nJ = length(cellChannelJrc); 
% jrc output: spikeTimeJrc, timeJrc, cluJrc, cellChannelJrc, wavJrc, nJ

%% kilosort-phy (you need https://github.com/kwikteam/npy-matlab)
spike_times = readNPY(fullfile(phy_dir, 'spike_times.npy')); % timestamp of all spikes, (n_spike, 1)
spike_template = readNPY(fullfile(phy_dir, 'spike_templates.npy')); % automated cluster of all spikes, (n_spike, 1)
spike_clusters = readNPY(fullfile(phy_dir, 'spike_clusters.npy'));  % final cluster of all spikes, (n_spike, 1)
template = readNPY(fullfile(phy_dir, 'templates.npy')); % template waveform of all clusters, (n_original_cluster, n_timepoint, all_valid_channel)

% I find the main channel of each template which has the lowest amplitude.
[~, mainSiteTemplate] = min(min(template, [], 2), [], 3);

% % load cluster data...
% cn_name = fullfile(phy_dir, 'cluster_groups.csv');
% fid = fopen(cn_name, 'r');
% cn = textscan(fid, '%f%s%[^\n\r]', 'Delimiter', '\t', 'TextType', 'string', 'Headerlines', 1, 'EndOfLine', '\r\n');
% fclose(fid);

% pick only good units
%inPhy = strcmp(cn{2}, 'good'); 
unitNumber = unique(spike_template)+1; 
%unitNumber = cn{1}(inPhy); % final good unit number list
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
%cn{1,1}(find(cellfun(@isempty,cluPhy)==1)) % there are some channels empty

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
for iC = 1:nC % increment source cells
    inSite = find(abs(cellChannelPhy - cellChannel(iC))<=4); % Phy clusters found near the current source cluster
    nIn = length(inSite);
    
    rhoPhy{iC} = zeros(nIn, 2);
    rhoPhy{iC}(:, 1) = inSite;
    
    iIndex = 1;
    for iI = 1:nIn % increment Phy clusters detected at nearby sites (<=4)
        rhoTemp = zeros(1, 30); % timelags 
        for iR = 1:30 % get corr between wavSrc and wavPhy
            if iR+51>size(wavPhy,2)
                tempWavPhy = interp1(iR:size(wavPhy,2),wavPhy(inSite(iI),iR:size(wavPhy,2)), linspace(iR,size(wavPhy,2),52)); 
                rhoTemp(iR) = abs(corr(wavSrc(10:61, iC), tempWavPhy')); % caution! the width of phy Wav equals 61! src and jrc Wavs are 63 (corresponding to 2.1ms).  
        
            end
        end
        rhoPhy{iC}(iI, 2) = max(rhoTemp);
        
        if rhoPhy{iC}(iI,2) >= 0.5 % if the phy cluster Wav is reasonably similar
            matchPhy{iC, iIndex} = zeros(nSpk(iC), 4, 'single');
            
            for iS = 1:nSpk(iC) % increment the current source spikes
                [minValue, minIndex] = min(abs(spikeTimePhy{inSite(iI)} - spikeTime{iC}(iS))); % spot the closest phy cluster spike 
                if minValue <= 15 % if there's a spike within this range (ms ?)
                    matchPhy{iC, iIndex}(iS, :) = [inSite(iI), minIndex, spikeTime{iC}(iS), spikeTimePhy{inSite(iI)}(minIndex)]; % phyClusterIdx, closestCurrPhySpkIdx, srcSpkTime, closestCurrPhySpkTime
                end
            end
            perfPhy(iC, iIndex) = nanmean(matchPhy{iC, iIndex}(:, 1)>0); % this basically tells how much in proportion the curr srcCluster (total 200) could be recovered by a phyCluster with a reasonably similar waveform. 
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
                if minValue <= 15 %15
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
for iC = 1:nC % increment the src cluster
    if ~isempty(matchPhy{iC, phyIndex(iC)}) % phyIndex indicates the best performing phy cluster among multiple possible phy clusters matched to the current srcCluster
        inCluPhy = find(matchPhy{iC, phyIndex(iC)}(:, 1)>0); % detected spikes, phyIndex: which matchPhyColumn
        if ~isempty(inCluPhy)
            targetClusterPhy = matchPhy{iC, phyIndex(iC)}(inCluPhy(1), 1); % which phy cluster
            nFound = sum(matchPhy{iC, phyIndex(iC)}(:, 1)>0); % phyIndex indicates the best performing phy cluster among multiple possible phy clusters matched to the current srcCluster
            falsePhy{targetClusterPhy} = [falsePhy{targetClusterPhy}, 1 - (nFound / nPhy(targetClusterPhy))];
        end
    end
    
    if ~isempty(matchJrc{iC, jrcIndex(iC)})
        inCluJrc = find(matchJrc{iC, jrcIndex(iC)}(:, 1)>0);
        if ~isempty(inCluJrc)
            targetClusterJrc = matchJrc{iC, jrcIndex(iC)}(inCluJrc(1), 1);
            nFound = sum(matchJrc{iC, jrcIndex(iC)}(:, 1)>0);
            falseJrc{targetClusterJrc} = [falseJrc{targetClusterJrc}, 1 - (nFound / nJrc(targetClusterJrc))];
        end
    end
end

falsePhy(cellfun(@isempty, falsePhy)) = {1}; 
falsePhyRate = 1-sum(cell2mat(cellfun(@(x) max(x < 0.1), falsePhy, 'UniformOutput', false)))/length(falsePhy);
falseJrc(cellfun(@isempty, falseJrc)) = {1};
falseJrcRate = 1-sum(cell2mat(cellfun(@(x) max(x < 0.1), falseJrc, 'UniformOutput', false)))/length(falseJrc);

subplot(1, 2, 1);
plot(-peakAmplitude, performancePhy, 'b*');
title('Kilosort2 performance'); xlabel('spikeAmp. uV'); ylim([0 1])

subplot(1, 2, 2);
plot(-peakAmplitude, performanceJrc, 'r*');
title('Jrc performance'); xlabel('spikeAmp. uV'); ylim([0 1])

subplot(1, 2, 1);
plot(frSrc, performancePhy, 'b*');
title('Kilosort2 performance'); xlabel('FR(Hz)'); ylim([0 1])

subplot(1, 2, 2);
plot(frSrc, performanceJrc, 'r*');
title('Jrc performance'); xlabel('FR (Hz)'); ylim([0 1])

subplot(1, 2, 1);
plot(-min(wavPhy, [], 2), cellfun(@min, falsePhy), 'b*'); % I don't know what's the phy's waveform unit is but it is very small.
title('Kilosort false positive / overmerged'); xlabel('spikeAmp. A.U.')

subplot(1, 2, 2);
plot(-min(squeeze(wavJrc(:, 1, :)), [], 1)', cellfun(@min, falseJrc), 'r*');
title('Jrc false positive / overmerged'); xlabel('spikeAmp. uV')

%% Analyze kilosort-dropped units 
missingCluPhy(:,1) = find(isnan(performancePhy) | performancePhy<0.7); % src clusters with poor spike reconstruction rate 
missingCluPhy(:,2) = cellChannel(missingCluPhy(:,1))'; % sites with poorly reconstructed clusters
missingCluPhy(:,3) = peakAmplitude(missingCluPhy(:,1));% peak amplitude of the Ks2 missing clusters
srtMissingCluPhy = sortrows(missingCluPhy,3); % 

spkTimeOffsetKs = cell(size(missingCluPhy,1),1);
for k = 1:size(missingCluPhy,1) % increment Ks2-missing clusters
    %plot(wavSrc(:,srtMissingCluPhy(k,1))) % missing src wav
    wavPhyCh = find(cellChannelPhy==srtMissingCluPhy(k,2)); % phy clu id recorded in that Ch (could be > 1)
    srcSpikeTimes = spikeTime{srtMissingCluPhy(k,1)}; % missing src wav spikeTimes
    
    for kk = 1:length(wavPhyCh)
        if ~isempty(wavPhyCh(kk)) && ~isempty(spikeTimePhy{wavPhyCh(kk)})
            %plot(wavPhy(wavPhyCh(kk),:))
            %length(spikeTimePhy{wavPhyCh(kk)})
            spkTimeOffsetKs{k,1} = zeros(length(srcSpikeTimes),2); % get the spike time offset
            for iS = 1:length(srcSpikeTimes) % increment the current source spikes
                [spkTimeOffsetKs{k,kk}(iS,1), spkTimeOffsetKs{k,kk}(iS,2)] = min(abs(spikeTimePhy{wavPhyCh(kk)} - srcSpikeTimes(iS))); % spot the closest phy cluster spike
            end
            
        end
    end
end
clearvars k




