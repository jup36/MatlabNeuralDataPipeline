filePaths = {
    '/Volumes/Extreme SSD/js2p0/WR38_052119/Matfiles', ... % Dual recording without Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR38_052419/Matfiles', ... % Corticostriatal recording M1 silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR39_100219/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR44_031020/Matfiles'};    % Dual recording with contra Cg delayed silencing (checked)

for f = 1:length(filePaths)
    jkvtPath = dir(fullfile(filePaths{f}, 'js2p0_tbytSpkHandJsTrjBin_50ms_stimPstimPrepExtWoTo_GPFA*')); 
    load(fullfile(jkvtPath.folder, jkvtPath.name), 'jkvt')

    trStartC = {jkvt.trStart}; 
    trJsReadyC = {jkvt.trJsReady}; 
    stimLaserOnC = {jkvt.stimLaserOn}; 
    stimLaserOffC = {jkvt.stimLaserOff}; 

    stimLaserOnRelToJsReadyC{f} = cell2mat(cellfun(@(a, b) a-b, stimLaserOnC, trJsReadyC, 'un', 0)); 
    stimLaserOnRelToTrStartC{f} = cell2mat(cellfun(@(a, b) a-b, stimLaserOnC, trStartC, 'un', 0)); 
    stimLaserDur{f} = cell2mat(cellfun(@(a, b) a-b, stimLaserOffC, stimLaserOnC, 'un', 0)); 

end

stimLaserOnRelToTrStartAllTrials = cell2mat(stimLaserOnRelToTrStartC); 
medianStimOnsetRelToTrStart = cellfun(@nanmedian, stimLaserOnRelToTrStartC); 
