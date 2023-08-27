function parse_tbyt_behavior(filePath)
% e.g., filePath = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA001/DA001_072623'; 

%% locate the nidq folder
filePath_nidq = GrabFiles_sort_trials('_g', 0, {filePath}); 
if isempty(filePath_nidq)
    filePath_nidq = uigetdir(filePath, 'Select a folder');
end

if length(filePath_nidq) > 1
    warning("More than 1 nidq folders found! Will run with the first one!")
end
filePath_nidq = filePath_nidq{1}; % just take the path

%% Load bin file
binFile = dir(fullfile(filePath_nidq, '*.nidq.bin')); % look for nidq.bin file
binFileI = cell2mat(cellfun(@(a) isletter(a(1)), {binFile.name}, 'un', 0));
binFile = binFile(binFileI);
if length(binFile)>1 || isempty(binFile)
    error('File could not be found or multiple nidq.bin files exist!');
end
binName = binFile.name;

%% Parse the corresponding metafile
meta  = ReadMeta(binName, filePath_nidq); % get the meta data (structure)
sample_rate=str2double(meta.niSampRate);
nSamp = floor(SampRate(meta));          % sampling rate 

%% Load evtInS or create one
evtInS = timestamp_behav_events(filePath_nidq, false, 'cmosExp', 'lick', 'faceCam', 'water', 'airpuff', 'photoDiode'); % behavioral events

%% tbytDat
numbTr = size(evtInS.photoDiode, 1); 
dat = struct; 
% stim On
stimOnC = num2cell(evtInS.photoDiode(:, 1)); 
[dat(1:numbTr).stimOn] = deal(stimOnC{:});  
% stim Off
stimOffC = num2cell(evtInS.photoDiode(:, 2)); 
[dat(1:numbTr).stimOff] = deal(stimOffC{:});  

% trial-by-trial peri-stim licks
tbytLickLogic = arrayfun(@(a) a-5 <= evtInS.lick & evtInS.lick <= a+5, [dat(:).stimOn], 'un', 0); 
tbytLicks = cellfun(@(a) evtInS.lick(a), tbytLickLogic, 'un', 0); 
[dat(1:numbTr).licks] = tbytLicks{:}; 

% trial-by-trial during-stim licks

% Create your 1x2 struct array
myStruct = struct('data', []);

% Use deal to fill the 'data' field with [1, 2]
[myStruct(1:2).data] = deal(1, 2);





end
