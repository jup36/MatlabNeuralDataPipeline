function individualUnitLaserEffectPlotGramm(filePath, dataS, unit)

unitIdFmt = 'Unit#%d'; 
% reach vs all reach-triggered laser trials
fig1 = spikeRasterGramm( [2e3 2e3], {'reach','stimLaser'}, [1.9e3 1.9e3], dataS.reach.SpkTimes{unit}, dataS.stmLaser.SpkTimes{unit} );
pbaspect([1 1 1])
%print(fig1, strcat(saveNameHeading,'reachVSlaser','Unit',num2str(unit)), '-dpdf'); % print figure as pdf

% pre-tagging vs tagging
fig2 = spikeRasterGramm( [5e3 5e3], {'taggingTrials'}, [2e3 2e3], dataS.tagLaser.SpkTimes{unit});
pbaspect([1 1 1])
%print(fig2, strcat(saveNameHeading,'tagVSpreTag','Unit',num2str(unit)), '-dpdf'); % print figure as pdf

if isfolder(fullfile(filePath,'Figure')) % if there's Figure folder already
    print(fig1, fullfile(filePath,'Figure',strcat('reachVSstimLaser',sprintf(unitIdFmt,unit))),'-dpdf','-bestfit')
    print(fig2, fullfile(filePath,'Figure',strcat('tagLaser',sprintf(unitIdFmt,unit))),'-dpdf','-bestfit')
else
    mkdir(fullfile(filePath,'Figure')) % otherwise, make a folder
    print(fig1, fullfile(filePath,'Figure',strcat('reachVSstimLaser',sprintf(unitIdFmt,unit))),'-dpdf','-bestfit')
    print(fig2, fullfile(filePath,'Figure',strcat('tagLaser',sprintf(unitIdFmt,unit))),'-dpdf','-bestfit')
end

end