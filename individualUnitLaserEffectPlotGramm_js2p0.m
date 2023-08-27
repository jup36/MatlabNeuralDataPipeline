function individualUnitLaserEffectPlotGramm_js2p0(filePath, tagLaser, unit)

unitIdFmt = 'Unit#%d'; 

% pre-tagging vs tagging
fig1 = spikeRasterGramm( [5e3 5e3], {'taggingTrials'}, [2e3 2e3], tagLaser.SpkTimes{unit});
pbaspect([1 1 1])
%print(fig1, strcat(saveNameHeading,'tagVSpreTag','Unit',num2str(unit)), '-dpdf'); % print figure as pdf

if isfolder(fullfile(filePath,'Figure')) % if there's Figure folder already
    print(fig1, fullfile(filePath,'Figure',strcat('tagLaser',sprintf(unitIdFmt,unit))),'-dpdf','-bestfit')
else
    mkdir(fullfile(filePath,'Figure')) % otherwise, make a folder
    print(fig1, fullfile(filePath,'Figure',strcat('tagLaser',sprintf(unitIdFmt,unit))),'-dpdf','-bestfit')
end

end