function saveCommonEventsImecNidq( filePath )%overwrite_imec, overwrite_bcs)
%SAVEEVENTRULESWITCH2 Reads IMEC binary data and makes event file
%   SAVEVENTRULESWITCH(BFL) reads BF file and makes event file using bin file. It also
%   extracts BCS behavior record and synchronize with IMEC sync data. All
%   units are saved in seconds.

%% default data location
%IMEC_DATA_PATH = 'H:\WR38_052119';
%USERNAME = getenv('USERNAME');
%BCS_DATA_PATH = 'H:\WR38_052119'; %['C:\Users\', USERNAME, '\OneDrive - Howard Hughes Medical Institute\src\vr\matlab\data\'];

%% Find bin data files
imec_file = fileSelectorExit(IMEC_DATA_PATH, 1, '*.ap.bin'); % fileSelectorExit doesn't prompt the user selection menu

%% Save event data
n_bin = length(imec_file);
for i_bin = 1:n_bin
    % data file name
    [bin_folder, ~] = fileparts(imec_file{i_bin});
    bin_folder_split = strsplit(bin_folder, '\');
    data_file = fullfile(bin_folder, [bin_folder_split{end}, '_data.mat']);
    fprintf('\n');
    disp(['======== ', imec_file{i_bin}, ' ========']);
    
    % find nidq data
    folder_join = strjoin(bin_folder_split(1:end-1), '\');
    nidq_file = dir([folder_join, '\', '*\*.nidq.bin']);
    
    %% run
    if ~(exist(data_file, 'file')==2)
        [event_imec, meta_imec] = readEventBin(imec_file{i_bin}); %save.readEventBin(imec_file{i_bin});
        saveImecEvent(data_file, event_imec, meta_imec);
        clear event_imec
        
        %saveBcsEvent(data_file, BCS_DATA_PATH);
        
        if ~isempty(nidq_file)
            [event_nidq, meta_nidq] = readEventBin(fullfile(nidq_file(1).folder, nidq_file(1).name));
            saveNidqEvent(data_file, event_nidq, meta_nidq);
            clear event_nidq
        else
            disp('No NIDQ data.');
        end
    else

    end
end
%slack('saveEvent done');

function saveImecEvent(data_file, event_data, meta)
clearvars trStart trEnd

% 0: trial start (note that matlab is 1-based)
% 7: trial end

trStartEvt = double(bitget(event_data, 1, 'uint16')); %[0; diff(double(bitget(event_data, 1, 'uint16')))]; 
[trStartIdx,~,~] = detecteventbythreshold(trStartEvt', meta.imSampRate, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'correctLongPulse', true); % trial Start
trStart.timeImec = trStartIdx./meta.imSampRate*1000; % time in msec

trEndEvt = double(bitget(event_data, 8, 'uint16'));   %[0; diff(double(bitget(event_data, 1, 'uint16')))]; 
[trEndIdx,~,~] = detecteventbythreshold(trEndEvt', meta.imSampRate, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'detectLater', trStartIdx(1), 'correctLongPulse', true); % trial Start
trEnd.timeImec = trEndIdx./meta.imSampRate*1000; % time in msec

disp(['Saving IMEC data to ', data_file]);
if exist(data_file, 'file')==2
    save(data_file, 'trStart', 'trEnd', '-append');
else
    save(data_file, 'trStart', 'trEnd');
end

function saveNidqEvent(data_file, event_data, meta)
%% load bcs mat file
load(data_file, 'trStart', 'trEnd');

%% process event_data
% 33: trial start
% 36: trial end
BYTE_PER_V = 2^15 / meta.niAiRangeMax;
BYTE_THRESHOLD_SYNC = 0.5 * BYTE_PER_V;

trStartEvt = double(event_data(:,33) >= BYTE_THRESHOLD_SYNC); %[0; diff(double(bitget(event_data, 1, 'uint16')))]; 
[trStartIdx,~,~] = detecteventbythreshold(trStartEvt', meta.niSampRate, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'correctLongPulse', true); % trial Start
trStart.timeNidq = trStartIdx./meta.niSampRate*1000; % time in msec

trEndEvt = double(event_data(:,36) >= BYTE_THRESHOLD_SYNC);   %[0; diff(double(bitget(event_data, 1, 'uint16')))]; 
[trEndIdx,~,~] = detecteventbythreshold(trEndEvt', meta.niSampRate, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'detectLater', trStartIdx(1), 'correctLongPulse', true); % trial Start
trEnd.timeNidq = trEndIdx./meta.niSampRate*1000; % time in msec

disp(['Saving NIDQ data to ', data_file]);
if exist(data_file, 'file')==2
    save(data_file, 'Sync', 'Trial', '-append');
else
    save(data_file, 'Sync', 'Trial');
end