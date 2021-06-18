%% Individual cell plot with gramm
filePath = '/Volumes/RAID2/parkj/NeuralData/ITphys/IT01_121317/Matfiles'; 
S=load('binSpkCountCTXIT01_121317.mat'); % load the corresponding binSpkCount structure
saveNameHeading = 'CTXIT01_121317'; 
stimE = load('binSpkCountCTXIT01_121317_stimE.mat'); 
stimE = stimE.('stimE'); 

unit = 37; % unit = 6, 10
% reach vs all reach-triggered laser trials
laserVSreach=stimE.meanLaser-stimE.meanReach;
laserVSreach(:,2)=1:length(laserVSreach);
laserVSreach=sortrows(laserVSreach,1);
fig1 = spikeRasterGrammSpikeTimes( [2e3 2e3], {'reach','stimLaser'}, [1.9e3 1.9e3], S.reach.SpkTimes{unit}, S.stmLaser.SpkTimes{unit} );
pbaspect([1 1 1])
print(fig1, fullfile(filePath,'Figure',strcat(saveNameHeading,'reachVSlaser','Unit',num2str(unit))),'-painters','-dpdf','-bestfit'); % print figure as pdf

% pre-tagging vs tagging
tagVSpreTag=stimE.sumTagStmOn-stimE.sumPreTagStmOn;
tagVSpreTag(:,2)=1:length(tagVSpreTag);
tagVSpreTag=sortrows(tagVSpreTag,1);
fig2 = spikeRasterGramm( [5e3 5e3], {'taggingTrials'}, [2e3 2e3], S.tagLaser.SpkTimes{unit});
pbaspect([1 1 1])
print(fig2, fullfile(filePath,'Figure',strcat(saveNameHeading,'tagVSpreTag','Unit',num2str(unit))),'-painters','-dpdf','-bestfit'); % print figure as pdf


spikeRasterGramm( [2e3 3e3], {'reach'}, [2e3 2.5e3], S.reach.SpkTimes{origUnit});
print(fullfile(filePath,'Figure',sprintf(indvUnitFigNameFmt,valUnit,origUnit)),'-dpdf','-bestfit')
individualUnitLaserEffectPlotGramm(filePath, S, origUnit) % individual unit plot depicting laser stim effect during session & tagging