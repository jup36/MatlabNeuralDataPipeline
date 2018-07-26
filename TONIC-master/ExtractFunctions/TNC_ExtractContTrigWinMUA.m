function [] = TNC_ExtractContTrigWinMUA()
% Show MUA for a given electrode using ContTrigWin data
% DEPENDENCY:
%     PRE-analyzed with TNC_ContTrigWins

numTrodes   = size(physStruct.phys.trode,2);
numTrials   = size(physStruct.phys.trode(1).spkData,1);
numSamps    = size(physStruct.phys.trode(1).beta,2);

currParams.smthParams.rise  = 1;
currParams.smthParams.decay = 20;
[currParams.kernel]         = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
physStruct.phys.currParams  = currParams;

for j=1:numTrodes
    for i=1:numTrials
        currTrialData   = physStruct.phys.trode(j).spkData(i,:);
        
        [events]                                            = TNC_EventDetect(currTrialData,30,3);
        [spikes]                                            = TNC_EventExtract(currTrialData,currTrialData,events.inds,[30,50]);
        spikes.inds                                         = events.inds;
        [aSpikes]                                           = TNC_EventAlign(spikes,10);
        
        physStruct.phys.trode(j).mua.raster.trial(i).inds   = events.inds;
        physStruct.phys.trode(j).mua.raster.trial(i).wfs    = aSpikes.Sh_wfs;
        physStruct.phys.trode(j).mua.raster.trial(i).xs     = aSpikes.Sh_x;
        
        delta                                               = zeros(1,numSamps);
        delta(round(events.inds./30))                       = 1;
        physStruct.phys.trode(j).mua.raster.delta(i,:)      = delta;
        physStruct.phys.trode(j).mua.raster.psth(i,:)       = conv(delta,currParams.kernel,'same');
        
        clear events spikes aSpikes
    end
end
