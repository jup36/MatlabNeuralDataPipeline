function [ csvData ] = getcsvfile
%This function locates generated csvData (col1: spike timing, col2: cluster id, col3: max site)
% and load it. 
%addpath(genpath('/Volumes/RAID2/parkj/jrclust_repo/JRCLUST')) % add jrclust-repo

fileList = dir('*.csv'); % identify the *.ap.bin file 

if length(fileList)==1 
    csvData = csvread(fileList.name);
else 
    error('Check the .csv file; one csv file must be present in the working folder!')
end

end

