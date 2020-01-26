function [ S_clu,viTime_spk,viSite_spk,viSite2_spk] = getjrcmatVar
%This function loads the 'S_Clu' (the structure containing spike
% information), and 'viTime_spk' (timestamp per spike - to be divided by the sampling rate for conversion to actual time).  
% 'viSite_spk' (site with the peak spike amplitude).
% in the '*jrc.mat' file.

fileList = dir('*_jrc.mat'); % identify the *.ap.bin file 

if length(fileList)==1 
    load(fileList.name,'S_clu','viTime_spk','viSite_spk','viSite2_spk');
else 
    error('Check the *jrc.mat file; one *_jrc.mat file must exist in the working folder!')
end

end

