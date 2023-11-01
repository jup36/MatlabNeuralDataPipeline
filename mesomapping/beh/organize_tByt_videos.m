function organize_tByt_videos(filePath)
%% get tbytDat
tbytDatPath = GrabFiles_sort_trials('tbytDat', 0, {fullfile(filePath, 'Matfiles')}); 
load(fullfile(tbytDatPath{1}), 'tbytDat')

%% get raw video files and csv bodypart coordinate files
% raw video files
filePath_raw_vid = findFolderWithString(filePath, '_vid');
if isempty(filePath_raw_vid)
    filePath_raw_vid = uigetdir(filePath); 
end
vidFiles = GrabFiles_sort_trials('behvid', 0, {filePath_raw_vid}); 
if numel(tbytDat) == numel(vidFiles)    
    vidTrialC = cellfun(@(a) regexp(a, '_([0-9.]+)\.', 'tokens'), vidFiles, 'un', 0); 
    vidTrial = cell2mat(cellfun(@(c) str2double(c{1}{1}), vidTrialC, 'un', 0));
    [~, vidI] = sort(vidTrial); 
    vidNameSorted = vidFiles(vidI); 
    [tbytDat(:).vidFile] = deal(vidNameSorted{:}); 
else
    warning('The number of raw video files does not match with the number of trials!')
end    
disp("Completed organizing the raw video files!")

% DLC bodypart coordinate csv files
filePath_cropped_vid = findFolderWithString(filePath, '_vid_cropped');
csvFiles = dir(fullfile(filePath_cropped_vid, 'behvid*.csv'));
% csv files
if numel(tbytDat) == numel(csvFiles)
    csvNames = {csvFiles.name}; 
    csvTimeTokenC = cellfun(@(a) regexp(a, '(\d{2}_\d{2}_\d{2})', 'tokens'), csvNames, 'un', 0); 
    csvTimeC  = cellfun(@(c) hmsToDateTime(c{1}{1}), csvTimeTokenC, 'un', 0); 
    [~, csvTimeI] = sort([csvTimeC{:}]); % to ensure the csv files are in correct timely order
    csvNamesSorted = csvNames(csvTimeI); 
    csvNamesSorted = cellfun(@(a) fullfile(filePath_cropped_vid, a), csvNamesSorted, 'un', 0); 
    [tbytDat(:).csvFile] = deal(csvNamesSorted{:}); 
else
    warning('The number of csv files does not match with the number of trials!') % If some videos are missing each videoes need to be assigned using the trial indices in their names (TO DO)
end
disp("Completed organizing the csv DLC coordinate files!")


% cropped video files
croppedVidFiles = dir(fullfile(filePath_cropped_vid, 'behvid*cropped.mp4'));
if numel(tbytDat) == numel(croppedVidFiles)
    cropVidNames = {croppedVidFiles.name}; 
    cropVidTimeTokenC = cellfun(@(a) regexp(a, '(\d{2}_\d{2}_\d{2})', 'tokens'), cropVidNames, 'un', 0); 
    cropVidTimeC  = cellfun(@(c) hmsToDateTime(c{1}{1}), cropVidTimeTokenC, 'un', 0); 
    [~, cropVidTimeI] = sort([cropVidTimeC{:}]); % to ensure the cropVid files are in correct timely order
    cropVidNamesSorted = cropVidNames(cropVidTimeI); 
    cropVidNamesSorted = cellfun(@(a) fullfile(filePath_cropped_vid, a), cropVidNamesSorted, 'un', 0); 
    [tbytDat(:).cropVidFile] = deal(cropVidNamesSorted{:}); 
else
    warning('The number of cropVid files does not match with the number of trials!') % If some videos are missing each videoes need to be assigned using the trial indices in their names (TO DO)
end
disp("Completed organizing the cropped video files!")

save(fullfile(tbytDatPath{1}), 'tbytDat')


end