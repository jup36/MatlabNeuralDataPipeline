function [unitTimeTrialBZ] = zscoreUnitTimeTrialB( base, unitTimeTrial )
% get the mean and std across all baseline time bins and trials, then z-score normalize unitTimeTrial 
base = nanmean(base,3); 
rsBase = reshape(base, [size(base,1), size(base,2)*size(base,3)]);
meanBase = nanmean(rsBase,2);
stdBase = nanstd(rsBase,0,2);

rsUnitTimeTrial = reshape(unitTimeTrial, [size(unitTimeTrial,1), size(unitTimeTrial,2)*size(unitTimeTrial,3)]); 
rsUnitTimeTrialBZ=(rsUnitTimeTrial-repmat(meanBase,[1 size(rsUnitTimeTrial,2)]))./repmat(stdBase,[1 size(rsUnitTimeTrial,2)]);
unitTimeTrialBZ = reshape(rsUnitTimeTrialBZ, [size(unitTimeTrial,1), size(unitTimeTrial,2), size(unitTimeTrial,3)]);

end