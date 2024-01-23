% 


filePath = '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles/js2p0_tbytSpkHandJsTrjBin_WR40_082019.mat'; 
load(fullfile(filePath), 'ss', 'jkvt')

stimTrI = cellfun(@(a) ~isempty(a), {ss.spkTimeBlaserI}); 

stbLaserIC = {ss(stimTrI).spkTimeBlaserI}; % spikeTimeBin laser index 
stbRchIC = {ss(stimTrI).spkRchIdx}; % spikeTimeBin reach index
stbPullIC = {ss(stimTrI).spkPullIdx}; % spikeTimeBin pull index


