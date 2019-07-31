function saveEventRuleSwitch2(overwrite_imec, overwrite_bcs)
%SAVEEVENTRULESWITCH2 Reads IMEC binary data and makes event file
%   SAVEVENTRULESWITCH(BFL) reads BF file and makes event file using bin file. It also
%   extracts BCS behavior record and synchronize with IMEC sync data. All
%   units are saved in seconds.

% overwrite: 1, yes; 0, no; 2, ask
if nargin < 1
    overwrite_imec = 0;
end
if nargin < 2
    overwrite_bcs = 0;
end

%% default data location
IMEC_DATA_PATH = 'E:';
USERNAME = getenv('USERNAME');
BCS_DATA_PATH = ['C:\Users\', USERNAME, '\OneDrive - Howard Hughes Medical Institute\src\vr\matlab\data\'];


%% Find bin data files
imec_file = fileSelector(IMEC_DATA_PATH, 1, '*.ap.bin');


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
    name_split = strsplit(bin_folder_split{end}, '_');
    name_join = strjoin(name_split(1:2), '_');
    folder_join = strjoin(bin_folder_split(1:end-1), '\');
    nidq_file = dir([folder_join, '\', name_join, '*\*.nidq.bin']);
    
    %% run
    if ~(exist(data_file, 'file')==2)
        [event_imec, meta_imec] = save.readEventBin(imec_file{i_bin});
        saveImecEvent(data_file, event_imec, meta_imec);
        clear event_imec
        
        saveBcsEvent(data_file, BCS_DATA_PATH);
        
        if ~isempty(nidq_file)
            [event_nidq, meta_nidq] = save.readEventBin(fullfile(nidq_file(1).folder, nidq_file(1).name));
            saveNidqEvent(data_file, event_nidq, meta_nidq);
            clear event_nidq
        else
            disp('No NIDQ data.');
        end
        
    else
        eventFileInfo = who('-file', data_file);
        
        checkVariable1 = {'Sync', 'Lick', 'Trial'};
        if ~all(ismember(checkVariable1, eventFileInfo)) || overwrite_imec == 1
            doCalcImec = 'Yes';
        elseif overwrite_imec == 2
            doCalcImec = questdlg('Imec event file already exists. Do you want to overwrite?', 'Overwrite', 'No');
        else
            doCalcImec = 'No';
        end
        
        if strcmp(doCalcImec, 'Yes')
            [event_imec, meta_imec] = save.readEventBin(imec_file{i_bin});
            saveImecEvent(data_file, event_imec, meta_imec);
            clear event_imec
        end
        
        checkVariable2 = {'Vr'};
        if ~all(ismember(checkVariable2, eventFileInfo)) || overwrite_bcs == 1
            doCalcBcs = 'Yes';
        elseif overwrite_bcs == 2
            doCalcBcs = questdlg('BCS event file already exists. Do you want to overwrite?', 'Overwrite', 'No');
        else
            doCalcBcs = 'No';
        end
        
        if strcmp(doCalcBcs, 'Yes')
            saveBcsEvent(data_file, BCS_DATA_PATH);
            
            if ~isempty(nidq_file)
                [event_nidq, meta_nidq] = save.readEventBin(fullfile(nidq_file(1).folder, nidq_file(1).name));
                saveNidqEvent(data_file, event_nidq, meta_nidq);
                clear event_nidq
            else
                disp('No NIDQ data.');
            end
        end
    end
end
slack('saveEvent done');




function saveImecEvent(data_file, event_data, meta)
clearvars Sync Trial Lick


% 0: sync
% 1: recording enable
% 2: trial
% 3: reward
% 4: noreward
% 5: x
% 6: x
% 7: lick

syncEvent = [0; diff(double(bitget(event_data, 1, 'uint16')))];
Sync.timeImec = find(syncEvent~=0)/meta.imSampRate;
Sync.typeImec = syncEvent(syncEvent~=0)==1;

lickEvent = [0; diff(double(bitget(event_data, 8, 'uint16')))];
Lick.timeImec = find(lickEvent~=0)/meta.imSampRate;
Lick.typeImec = lickEvent(lickEvent~=0)==1;

trialEvent = [0; diff(double(bitget(event_data, 3, 'uint16')))];
Trial.timeStartImec = find(trialEvent == 1)/meta.imSampRate;
Trial.timeDelayImec = find(trialEvent == -1)/meta.imSampRate;

trialCorrectEvent = [0; diff(double(bitget(event_data, 4, 'uint16')))];
trialWrongEvent = [0; diff(double(bitget(event_data, 5, 'uint16')))];
Trial.timeCorrectImec = find(trialCorrectEvent==1)/meta.imSampRate;
Trial.timeWrongImec = find(trialWrongEvent==1)/meta.imSampRate;

Trial.nCorrectImec = length(Trial.timeCorrectImec);
Trial.nWrongImec = length(Trial.timeWrongImec);
Trial.nTrialImec = Trial.nCorrectImec + Trial.nWrongImec;

timeResult = [Trial.timeCorrectImec; Trial.timeWrongImec];
[Trial.timeResultImec, trialIndex] = sort(timeResult);
Trial.resultImec = trialIndex <= Trial.nCorrectImec;

disp(['Saving IMEC data to ', data_file]);
if exist(data_file, 'file')==2
    save(data_file, 'Sync', 'Lick', 'Trial', '-append');
else
    save(data_file, 'Sync', 'Lick', 'Trial');
end




function data_file = saveBcsEvent(data_file, bcs_path)
%% find bcs file
[file_folder, file_name] = fileparts(data_file);
file_name_split = strsplit(file_name, '_');

bin_file = dir(fullfile(file_folder, '*.imec.ap.bin'));
bcs_file = dir(fullfile(bcs_path, [file_name_split{1},'_',file_name_split{2},'_*.mat']));

n_bcs_file = length(bcs_file);
time_diff = zeros(n_bcs_file, 1);
for i_file = 1:n_bcs_file
    time_diff(i_file) = abs(bcs_file(i_file).datenum - bin_file.datenum);
end
[min_time_diff, index_time_diff] = min(time_diff);
if min_time_diff < 2/24/60 % 2 minite diff
    session_path = bcs_path;
    session_file = bcs_file(index_time_diff).name;
else
    [session_file, session_path] = uigetfile([bcs_path, file_name_split{1},'_',file_name_split{2},'_*.mat'], 'Choose session file');
    if session_file==0; return; end
end


%% load bcs mat file
load(fullfile(session_path, session_file), 'data', 'nTrial');
load(data_file, 'Sync', 'Trial', 'Lick');


% trial number
Trial.nTrialBcs = nTrial;
if Trial.nTrialBcs ~= Trial.nTrialImec
    error('Trial number mismatch!');
end

%% sync check
nSyncImec = length(Sync.typeImec);
nSyncBcs = length(data.syncType);
disp(['IMEC sync pulse: ', num2str(nSyncImec),', BCS sync pulse: ', num2str(nSyncBcs)]);

assert(nSyncImec - nSyncBcs <= 1, 'Sync number mismatch');
nSync = min(nSyncImec, nSyncBcs);
assert(all(Sync.typeImec(1:nSync)==data.syncType(1:nSync)), 'Disrupted Sync!');
Sync.timeImec = Sync.timeImec(1:nSync);
Sync.typeImec = Sync.typeImec(1:nSync);
Sync.timeBcs = data.syncTime(1:nSync);
Sync.typeBcs = data.syncType(1:nSync);


% sync, check time overflow
Sync.timeBcs = double(Sync.timeBcs);
subplot(2, 2, 1); plot(Sync.timeBcs);
syncOverflow = find(diff(Sync.timeBcs) < -4E9) + 1;
if ~isempty(syncOverflow)
    disp(['Sync overflow at ', num2str(syncOverflow')]);
    for iS = 1:length(syncOverflow)
        Sync.timeBcs(syncOverflow(iS):end) = Sync.timeBcs(syncOverflow(iS):end) + double(intmax('uint32'));
    end
end
Sync.timeBcs = Sync.timeBcs / 10^6;
subplot(2, 2, 2); plot(Sync.timeBcs);

% Vr, check time overflow
Vr.timeBcs = double(data.timeStamp(:, 1));
subplot(2, 2, 3); plot(Vr.timeBcs);
timeOverflow = find(diff(Vr.timeBcs) < -4E9) + 1;
if ~isempty(timeOverflow)
    disp(['Time overflow at ', num2str(timeOverflow')]);
    for iT = 1:length(timeOverflow)
        Vr.timeBcs(timeOverflow(iT):end) = Vr.timeBcs(timeOverflow(iT):end) + double(intmax('uint32'));
    end
end
Vr.timeBcs = Vr.timeBcs / 10^6;
subplot(2, 2, 4); plot(Vr.timeBcs);


% Lick time overflow
Lick.timeBcs = double(data.lickTime);
lickOverflow = find(diff(Lick.timeBcs) < -4E9) + 1;
if ~isempty(lickOverflow)
    disp(['Lick overflow at ', num2str(lickOverflow')]);
    for iT = 1:length(lickOverflow)
        Lick.timeBcs(lickOverflow(iT):end) = Lick.timeBcs(lickOverflow(iT):end) + double(intmax('uint32'));
    end
end
Lick.timeBcs = Lick.timeBcs / 10^6;
Lick.typeBcs = data.lickType == 1;


%% Vr time, data
Vr.timeImec = interp1(Sync.timeBcs, Sync.timeImec, Vr.timeBcs, 'linear', 'extrap');
Vr.position = double(data.position(:,1:2))/100; % in centimeter
Vr.speed = double(data.velocity)/100; % cm/s
Vr.ballVelocity = data.ballvelocity; % cm/s
Vr.roll = data.roll;
Vr.pitch = data.pitch;
% Vr.yaw = data.yaw;
Vr.event = data.event;


%% Trial time
trialTime = double(data.trialTime);
trialOverflow = find(diff(trialTime) < -2E9) + 1;
if ~isempty(trialOverflow)
    disp(['Trial time overflow at ', num2str(trialOverflow')]);
    for iT = 1:length(trialOverflow)
        trialTime(trialOverflow(iT):end) = trialTime(trialOverflow(iT):end) + double(intmax('uint32'));
    end
end
trialTime = trialTime / 10^6;

inTrial = data.trial(:, 1)==1;
timeStart = trialTime(inTrial);
Trial.timeStartBcs = interp1(Sync.timeBcs, Sync.timeImec, timeStart(1:Trial.nTrialBcs), 'linear', 'extrap');


% Teleportation time
teleportIndex = [false; (abs(diff(Vr.position(:, 1))) >= 40 | abs(diff(Vr.position(:, 2))) >= 100)];
nTimeStart = length(Trial.timeStartBcs);
Trial.timeStartVr = NaN(nTimeStart, 1);

for iT = 1:nTimeStart
    inTrialTime = Vr.timeImec >= Trial.timeStartBcs(iT) & Vr.timeImec < Trial.timeStartBcs(iT) + 1;
    timeStartVrIndex = find(inTrialTime & teleportIndex, 1, 'first');
    if isempty(timeStartVrIndex); continue; end
    Trial.timeStartVr(iT) = Vr.timeImec(timeStartVrIndex);
    if Trial.timeStartVr(iT) > Trial.timeStartBcs(iT) + 0.2
        disp(['Delayed teleport at trial ', num2str(iT)]);
    end
end

nTimeDelay = length(Trial.timeResultImec);
Trial.timeDelayVr = NaN(nTimeDelay, 1);

for iT = 1:nTimeDelay
    inTrialTime = Vr.timeImec >= Trial.timeResultImec(iT) & Vr.timeImec < Trial.timeResultImec(iT) + 1;
    timeDelayVrIndex = find(inTrialTime & teleportIndex, 1, 'first');
    if isempty(timeDelayVrIndex); continue; end
    Trial.timeDelayVr(iT) = Vr.timeImec(timeDelayVrIndex);
end

% Trial cue
inTrial = data.trial(:, 1) == 2 | data.trial(:, 1) == 3;
Trial.resultBcs = data.trial(inTrial, 1) == 2;

assert(all(Trial.resultBcs == Trial.resultImec));

cueType = data.trial(inTrial, 3);
timeResult = trialTime(inTrial);
Trial.timeResultBcs = interp1(Sync.timeBcs, Sync.timeImec, timeResult, 'linear', 'extrap');

Trial.task = floor(cueType / 100);

[Trial.cue, Trial.target, Trial.choice] = deal(NaN(Trial.nTrialBcs, 1));
for iT = 1:Trial.nTrialBcs
    switch cueType(iT)            
        case 101
            Trial.cue(iT) = NaN;
            Trial.target(iT) = 1;
            
        case 102
            Trial.cue(iT) = NaN;
            Trial.target(iT) = 2;
            
        case 121
            Trial.cue(iT) = 1;
            Trial.target(iT) = 1;
            
        case 122
            Trial.cue(iT) = 2;
            Trial.target(iT) = 2;
            
        case 123
            Trial.cue(iT) = 2;
            Trial.target(iT) = 1;
            
        case 124
            Trial.cue(iT) = 1;
            Trial.target(iT) = 2;
            
        case 201
            Trial.cue(iT) = 1;
            Trial.target(iT) = 1;
            
        case 202
            Trial.cue(iT) = 2;
            Trial.target(iT) = 2;
    end
    
    if Trial.resultBcs(iT)
        Trial.choice(iT) = Trial.target(iT);
    else
        Trial.choice(iT) = 3 - Trial.target(iT);
    end
    
    if cueType(iT) == 100
        Trial.cue(iT) = NaN;
        Trial.target(iT) = NaN;
        
        leftChoiceIndex = strncmp(data.event, 'l', 1);
        rightChoiceIndex = strncmp(data.event, 'r', 1);
        choiceIndex = leftChoiceIndex | rightChoiceIndex;
        inTimeResult = (Vr.timeImec > Trial.timeResultBcs(iT) - 0.2 & Vr.timeImec <= Trial.timeResultBcs(iT));
        
        Trial.choice(iT) = rightChoiceIndex(find(inTimeResult & choiceIndex, 1, 'first')) + 1;
    end
end
            
Trial.pChoice = [NaN; Trial.choice(1:end-1)];


disp(['Saving BCS data to ', data_file]);
if exist(data_file, 'file')==2
    save(data_file, 'Sync', 'Vr', 'Trial', 'Lick', '-append');
else
    save(data_file, 'Sync', 'Vr', 'Trial', 'Lick');
end



function saveNidqEvent(data_file, event_data, meta)
%% load bcs mat file
load(data_file, 'Sync', 'Trial');

%% process event_data
% 1: photodiode data
% 2: sync data
BYTE_PER_V = 2^15 / meta.niAiRangeMax;
BYTE_THRESHOLD_SYNC = 0.5 * BYTE_PER_V;

syncEvent = [0; diff(event_data(:, 2) >= BYTE_THRESHOLD_SYNC)];
Sync.timeNidq = find(syncEvent~=0)/meta.niSampRate;
Sync.typeNidq = syncEvent(syncEvent~=0)==1;

% check if the number of sync pulses are the same
fprintf('IMEC sync: %d, BCS: %d, NIDQ: %d\n', length(Sync.timeImec), length(Sync.timeBcs), length(Sync.timeNidq));
assert(length(Sync.timeNidq) - length(Sync.timeImec) <= 1, 'Lengths of sync pulses are not the same.')

n_pulse = min([length(Sync.timeNidq), length(Sync.timeImec)]);
Sync.timeNidq = Sync.timeNidq(1:n_pulse);
Sync.typeNidq = Sync.typeNidq(1:n_pulse);

assert(all(Sync.typeNidq == Sync.typeImec), 'Pulse sign is wrong.');

plot(Sync.timeImec, Sync.timeNidq);

% find the exact trial start point
index_start_nidq = zeros(Trial.nTrialImec, 1);
sample_start_nidq = floor(interp1(Sync.timeImec, Sync.timeNidq, Trial.timeStartBcs, 'linear', 'extrap') * meta.niSampRate);
for i_trial = 1:Trial.nTrialImec
    base_zone = sample_start_nidq(i_trial) + (0:ceil(meta.niSampRate * 0.05));
    vr_threshold = max(event_data(base_zone, 1)) * 1.2;
    test_zone = sample_start_nidq(i_trial) + (floor(meta.niSampRate * 0.05):ceil(meta.niSampRate * 0.25));
    
    index_start_nidq(i_trial) = test_zone(find(event_data(test_zone, 1) > vr_threshold, 1, 'first'));
end

Trial.timeStartNidq = interp1(Sync.timeNidq, Sync.timeImec, index_start_nidq / meta.niSampRate, 'linear', 'extrap');
% The interquartile range were 0.1628 and 0.1810 seconds.

disp(['Saving NIDQ data to ', data_file]);
if exist(data_file, 'file')==2
    save(data_file, 'Sync', 'Trial', '-append');
else
    save(data_file, 'Sync', 'Trial');
end