function [trials] = getSuccessfulStimTrialsWithVideo(filePath)
    cd(fullfile(filePath))
    load(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles.mat'),'jsTime1k_KV')
    successStimTrials = find(~isnan([jsTime1k_KV(:).stimLaserOn]) & ~isnan([jsTime1k_KV(:).stimLaserOff]) & cellfun(@(a) strcmpi(a,'sp'), {jsTime1k_KV(:).trialType})==1); 
    
    
   

end