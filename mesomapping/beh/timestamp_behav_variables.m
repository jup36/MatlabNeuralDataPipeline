function timestamp_behav_events(path_raw, varargin)
%This function takes the number of behavioral variables as varargin and
% gets the timestamps (in sec) of the variables and saves them as a
% structure named evtInS in the files origianl nidq directory. 
% TO DO: Add more variables in the SWITCH CASE block. 

% load bin file
binFile = dir(fullfile(path_raw, '*.nidq.bin')); % look for nidq.bin file
binFileI = cell2mat(cellfun(@(a) isletter(a(1)), {binFile.name}, 'un', 0));
binFile = binFile(binFileI);

if length(binFile)>1 || isempty(binFile)
    error('File could not be found or multiple nidq.bin files exist!');
end

binName = binFile.name;

% Parse the corresponding metafile
meta  = ReadMeta(binName, path_raw); % get the meta data (structure)
sample_rate=str2double(meta.niSampRate);
nSamp = floor(SampRate(meta));          % sampling rate (default: 25kHz)
%totalTimeSecs = str2double(meta.fileTimeSecs); % total duration of file in seconds
%channels = textscan(meta.acqMnMaXaDw,'%n %n %n %n','Delimiter',',');

output = struct; 
evtInS = struct; 

% Main (load raw traces and perform event detection)
if exist(fullfile(path_raw,'gainCorrectRawTraces.mat'), 'file')==2
    for jj = 1:length(varargin)
        fieldName = sprintf('variable_%d', jj); 
        output.(fieldName) = load(fullfile(path_raw,'gainCorrectRawTraces.mat'), varargin{jj}); 
        assert(isfield(output.(fieldName), varargin{jj})) % ensure the raw trace is loaded properly
       
        switch varargin{jj}
            case 'lick'
                lick = output.(fieldName).lick; 
                timeStamp = 1:length(lick);
                timeStampSec = timeStamp./sample_rate;
                lickIdx = detecteventbythreshold(lick, nSamp, 50, 'stdFactor', 3, 'plotRez', true, 'chunkPulses', false, 'correctLongPulse', true); 
                evtInS.lick = timeStampSec(lickIdx); 
                fprintf('completed lick detection!\n');
            case 'faceCam'
                faceCam = output.(fieldName).faceCam; 
                timeStamp = 1:length(faceCam);
                timeStampSec = timeStamp./sample_rate;
                [faceCamRiseIdx, ~, faceCamTrainIdx] = detecteventbythreshold_noAbs(faceCam, nSamp, 2, 'stdFactor', 2.5, 'plotRez', true, 'chunkPulses', true, 'chunkInterval', 2000, 'correctLongPulse', false); % camera trigger
                evtInS.faceCam = [timeStampSec(faceCamRiseIdx)', faceCamTrainIdx']; 
                fprintf('completed faceCam detection!\n');
        end    
    end
end

%% save evtInS
if exist(fullfile(path_raw, "evtInS"), 'file') == 2
    save(fullfile(path_raw, "evtInS"), 'evtInS', '-append')
else
    save(fullfile(path_raw, "evtInS"), "evtInS")
end
fprintf('completed event detection and saved events in the nidq folder!\n');

end