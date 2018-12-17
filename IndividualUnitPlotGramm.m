%% Individual cell plot with gramm
S=load('binSpkCountCTXPT08_062218.mat'); % load the corresponding binSpkCount structure
saveNameHeading = 'CTXPT08_062218'; 

unit = 20;
% reach vs all reach-triggered laser trials
laserVSreach=stimE.meanLaser-stimE.meanReach;
laserVSreach(:,2)=1:length(laserVSreach);
laserVSreach=sortrows(laserVSreach,1);
fig1 = spikeRasterGramm( [2e3 2e3], {'reach','stimLaser'}, [1.9e3 1.9e3], S.reach.SpkTimes{unit}, S.stmLaser.SpkTimes{unit} );
pbaspect([1 1 1])
print(fig1, strcat(saveNameHeading,'reachVSlaser','Unit',num2str(unit)), '-dpdf'); % print figure as pdf

% pre-tagging vs tagging
tagVSpreTag=stimE.sumTagStmOn-stimE.sumPreTagStmOn;
tagVSpreTag(:,2)=1:length(tagVSpreTag);
tagVSpreTag=sortrows(tagVSpreTag,1);
fig2 = spikeRasterGramm( [5e3 5e3], {'taggingTrials'}, [2e3 2e3], S.tagLaser.SpkTimes{unit});
pbaspect([1 1 1])
print(fig2, strcat(saveNameHeading,'tagVSpreTag','Unit',num2str(unit)), '-dpdf'); % print figure as pdf
