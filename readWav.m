% Opens raw waveform from *_spkwav.jrc file.

% _spkwav.jrc file is a binary file containing filtered waveform per spike.
% Dimension is described in dimm_spk of _jrc.mat file
% Format: nSamples * nSites_spk * nSpikes = 32 * 14 * nSpike

% _spkraw.jrc file is a file containing raw waveforms per spike.
% Dimension is described in dimm_raw of _jrc.mat file
% Format: nSamples * nSites_spk * nSpikes

% Select file

IMEC_DATA_PATH = 'F:\neozig\data_ephys\';

dataFileList = dir(fullfile(IMEC_DATA_PATH, '*_spkwav.jrc'));
nData = length(dataFileList);
dataFile = cell(nData, 1);
for iData = 1:nData
    dataFile{iData} = fullfile(dataFileList(iData).folder, dataFileList(iData).name);
end

selection = listdlg('PromptString', 'Select a file', ...
    'SelectionMode', 'single', ...
    'ListSize', [550, 400], ...
    'ListString', dataFile); % create list selection dialog box

if isempty(selection); return; end
dataSelected = dataFile{selection};

% File size
jrcFile = replace(dataSelected, 'spkwav.jrc', 'jrc.mat'); 
load(jrcFile, 'dimm_spk', 'S_clu'); % load the 'jrc.mat' file

dataSize = dataFileList(selection).bytes;
nSpike = dataSize / 2 / 32 / 14; % nSamples * nSites_spk * nSpikes = 32 * 14 * nSpike

if dimm_spk(3) ~= nSpike
    disp('Something is wrong...');
    return;
end

% Open file
tic;
fid = fopen(dataSelected, 'rb');
dataArray = fread(fid, [32*14, nSpike], 'int16');
fclose(fid);
toc;

% Plot spikes
spikeIndex = S_clu.viClu;

iUnit = 1;
spike = dataArray(1:32, spikeIndex==iUnit); % just show main channel

plot(spike(:, 1:100));
