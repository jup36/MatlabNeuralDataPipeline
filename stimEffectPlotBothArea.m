function stimEffectPlotBothArea(filePath, binSpkFileName, stimEctx, stimEstr)
%This function is to quantify the effect of laser stimulation on neural
% activity of both the STR and CTX. 
% 'binSpkFileName' refers to the binSpkCount file containing binned PSTHs
% (e.g. 'binSpkCountStrCtxIT01_121317.mat'). 

cd(filePath)

% load the binSpkCount files (PSTHs)
binSpkFile = dir(fullfile(filePath, binSpkFileName));
S = load(binSpkFile.name); 

% load the 'binSpkCount*_stimE.mat'
stimEfiles = dir(fullfile(filePath,'*stimE.mat')); 

if length(stimEfiles)>2 % usually there are two '*_stimE.mat' files; CTX and STR 
    return
end
   

for f = 1:length(stimEfiles)
    
end


%% combine str & ctx variables 
meanLaser = [stimEstr.meanLaser; stimEctx.meanLaser];   % mean laser activity
meanReach = [stimEstr.meanReach; stimEctx.meanReach ];  % mean reach activity
meanTagStmOn = [stimEstr.meanTagStmOn; stimEctx.meanTagStmOn];          % mean Tag-stim activity 
meanPreTagStmOn = [stimEstr.meanPreTagStmOn; stimEctx.meanPreTagStmOn]; % mean pre Tag-stim activity

sites = cell2mat(S.reach.Site); % sites on the IMEC probe
FRidx = [stimEstr.FRidx; stimEctx.FRidx]; % firing rate index 

% plot meanLaser-meanReach
laserVSreach=meanLaser-meanReach;
laserVSreach(:,2)=1:length(laserVSreach);
laserVSreach=sortrows(laserVSreach,1);

fig1 = imecOpt3GeomColorMapZ( sites(laserVSreach(FRidx,2)), laserVSreach(FRidx,1), 'rb', true, stimEctx.p.Results.colorAxis, stimEctx.p.Results.probeDepth ); % imecOpt3GeomColorMapZ( imecSites, Z, colorScheme, drawAllSites, colorAxis, varargin )
print(fig1, strcat('binSpkCountStrCtx_stimE_siteColorMap'), '-dpdf'); % print figure as pdf

% plot meanTagStmOn-meanPreTagStmOn 
tagVSpreTag = meanTagStmOn - meanPreTagStmOn; 
tagVSpreTag(:,2)=1:length(tagVSpreTag);
tagVSpreTag =sortrows(tagVSpreTag,1);

fig2 = imecOpt3GeomColorMapZ( sites(tagVSpreTag(FRidx,2)), tagVSpreTag(FRidx,1), 'rb', true, stimEctx.p.Results.colorAxis, stimEctx.p.Results.probeDepth ); % imecOpt3GeomColorMapZ( imecSites, Z, colorScheme, drawAllSites, colorAxis, varargin )
print(fig2, strcat('binSpkCountStrCtx_tagE_siteColorMap'), '-dpdf'); % print figure as pdf

%% individual unit plot (better use the combined 'binSpkCountStrCtx*'.mat)
unit = 216;
fig1 = spikeRasterGramm( [1e3 3e3], {'reach','stimReach'}, [1e3 2e3], S.reach.SpkTimes{unit}, S.stmReach.SpkTimes{unit} );
print(fig1, strcat('reachVSstimReach','Unit',num2str(unit)), '-dpdf'); % print figure as pdf

% fig2 = spikeRasterGramm( [1e3 3e3], {'reach','stimLaser'}, [1e3 2e3], S.reach.SpkTimes{unit}, S.stmLaser.SpkTimes{unit} );
% print(fig2, strcat('reachVSlaser','Unit',num2str(unit)), '-dpdf'); % print figure as pdf

% pre-tagging vs tagging
fig3 = spikeRasterGramm( [5e3 5e3], {'taggingTrials'}, [1e3 1e3], S.tagLaser.SpkTimes{unit});
print(fig3, strcat('preTagVsTag','Unit',num2str(unit)), '-dpdf'); % print figure as pdf



end

