function [trStartImec, trEndImec] = saveImecEvent(data_file, event_data, meta)
clearvars trStart trEnd

% 0: trial start (note that matlab is 1-based)
% 7: trial end

trStartEvt = double(bitget(event_data, 1, 'uint16')); %[0; diff(double(bitget(event_data, 1, 'uint16')))]; 
[trStartIdx,~,~] = detecteventbythreshold(trStartEvt', meta.imSampRate, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'correctLongPulse', true); % trial Start
trStartImec = trStartIdx./meta.imSampRate*1000; % time in msec

trEndEvt = double(bitget(event_data, 8, 'uint16'));   %[0; diff(double(bitget(event_data, 1, 'uint16')))]; 
[trEndIdx,~,~] = detecteventbythreshold(trEndEvt', meta.imSampRate, 50, 'stdFactor', 1, 'plotRez', false, 'chunkPulses', false, 'detectLater', trStartIdx(1), 'correctLongPulse', true); % trial Start
trEndImec = trEndIdx./meta.imSampRate*1000; % time in msec

disp(['Saving IMEC data to ', data_file]);
if exist(data_file, 'file')==2
    save(data_file, 'trStartImec', 'trEndImec', '-append');
else
    save(data_file, 'trStartImec', 'trEndImec');
end
end