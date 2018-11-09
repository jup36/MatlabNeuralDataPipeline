function [fileNameMap] = renameAviFiles(fileDir, animalId, varargin)
%This function is to rename Avi files to be suited for
% 'setUpDirUsingCluster'.

tbytFileList = dir(fullfile(fileDir));    % trial-by-trial files
tbytFileList([tbytFileList.isdir]) = [];  % just take files
fileNames = {tbytFileList.name};

formatSpecF = 'M%03d_%8d_front_v%03d';
formatSpecS = 'M%03d_%8d_side_v%03d';
dateFormatOut = 'yyyymmdd'; 

for i = 1:length(fileNames)  %# Loop over the file names
    dateNum = str2double(datestr(tbytFileList(i).date,dateFormatOut)); 
    [~,~,ext] = fileparts(fullfile(fileDir,fileNames{i}));
    if ~isempty(strfind(tbytFileList(i).name, 'cam0'))
        newFileName = strcat(sprintf(formatSpecF, animalId, dateNum, i),ext);
    elseif ~isempty(strfind(tbytFileList(i).name, 'cam1')) 
        newFileName = strcat(sprintf(formatSpecS, animalId, dateNum, i),ext);
    else
    end        
    tbytFileList(i).nFileName = newFileName;
    f = fullfile(fileDir, newFileName);
    g = fullfile(fileDir, fileNames{i});
    movefile(g,f);        %# Rename the file
end

save( , tbytFileList)

