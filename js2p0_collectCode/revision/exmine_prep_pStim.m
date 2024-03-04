

reachTrialI = ~cell2mat(cellfun(@isempty, {jkvt.rStartToPull}, 'UniformOutput', false))'; 
stimTrialI = ~cell2mat(cellfun(@isnan, {jkvt.stimLaserOn}, 'UniformOutput', false))'; 

prepAlignC = cellfun(@(a) a-3000, {jkvt.rStartToPull}, 'UniformOutput', false);
prepAlignCIsempty = cell2mat(cellfun(@isempty, prepAlignC, 'UniformOutput', false)); 
prepAlignC(prepAlignCIsempty) = deal({NaN});
prepAlignC(stimTrialI) = deal({NaN}); 
prepAlign = cell2mat(cellfun(@(a, b) a-b, prepAlignC, {jkvt.trJsReady}, 'UniformOutput', false))'; 

pStimAlignC = {jkvt.pLaserOn};
pStimAlignCIsempty = cell2mat(cellfun(@isempty, pStimAlignC, 'UniformOutput', false)); 
pStimAlignC(pStimAlignCIsempty) = deal({NaN}); 
pStimAlignC(~reachTrialI) = deal({NaN}); 
pStimAlign = cell2mat(cellfun(@(a, b) a-b, pStimAlignC, {jkvt.trJsReady}, 'UniformOutput', false))'; 

prepPstimAlign = [prepAlign, pStimAlign]; 

missingPstimTrials = find(~isnan(prepAlign) & isnan(pStimAlign)); % this mismatch is due to 


rStartToPullC = {jkvt.rStartToPull}; 
rStartToPullCempty = cellfun(@isempty, rStartToPullC); 
rStartToPullC(rStartToPullCempty) = deal({NaN}); 

cellfun(@(a, b) a-b, rStartToPullC, {jkvt.trJsReady}, 'UniformOutput', false); 