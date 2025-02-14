function targdimred_extend_KN(fileDir, targVarName)
cd(fileDir)
targDir = dir(fullfile(fileDir,targVarName));
load(fullfile(targDir.folder, targDir.name), 'projectMat', 'unitTimeTrialBZ', 'unitTimeTrialStmLaser', 'unitIdx', 'ts', 'trialId', 'gaussianKernel2')

%% project to the new projection matrix
for t = 1:size(unitTimeTrialBZ,3) % trial-by-trial
    uttBZ = unitTimeTrialBZ(unitIdx,:,:); % unit indexing
    tmpNtj = projectMat'*uttBZ(:,:,t); % projection
    nTjCellAllTrials{1,1,t} = cell2mat(arrayfun(@(ROWIDX) conv(tmpNtj(ROWIDX,:),gaussianKernel2,'same'), (1:size(tmpNtj,1))', 'UniformOutput', false)); % smooth
    for uu = 1:size(uttBZ,1)
        tmpNtjEach = projectMat(uu,:)'*uttBZ(uu,:,t); % projection of each neuron
        nTjEachAllTrials{uu,1,t} = cell2mat(arrayfun(@(ROWIDX) conv(tmpNtjEach(ROWIDX,:),gaussianKernel2,'same'), (1:size(tmpNtjEach,1))', 'UniformOutput', false)); % smooth
    end
end
clearvars t
% separate out trajectories of perturb vs. non-perturb trials
valStimTrI = ts.reachNoStimIdx(trialId);
nTjC.nPtb = nanmean(reshape([nTjCellAllTrials{valStimTrI}],3,100,[]),3);
nTjC.ptb = nanmean(reshape([nTjCellAllTrials{~valStimTrI}],3,100,[]),3);

% trajectories of perturb vs. non-perturb trials per unit
for uu = 1:size(nTjEachAllTrials,1)
    prjKN_unit.nPtb{uu,1} = cell2mat(nTjEachAllTrials(uu,:,valStimTrI)); 
    prjKN_unit.ptb{uu,1} = cell2mat(nTjEachAllTrials(uu,:,~valStimTrI)); 
end

% trial-by-trial projection scores
for dd = 1:3
    prjKN_nPtb{dd,1} = cell2mat(squeeze(cellfun(@(a) a(dd,:), nTjEachAllTrials(valStimTrI), 'un', 0)));
    prjKN_ptb{dd,1} = cell2mat(squeeze(cellfun(@(a) a(dd,:), nTjEachAllTrials(~valStimTrI), 'un', 0)));
end

%% project the activity during stmLaser onto the new projection matrix
unitTimeTrialB_stmLaser = bin1msSpkCountMat( unitTimeTrialStmLaser, 50, 50, 'align', 'center' );
rsUnitTimeTrial_stmLaser = reshape(unitTimeTrialB_stmLaser, [size(unitTimeTrialB_stmLaser,1), size(unitTimeTrialB_stmLaser,2)*size(unitTimeTrialB_stmLaser,3)]);
meanPerUnit = nanmean(rsUnitTimeTrial_stmLaser,2);
stdPerUnit = nanstd(rsUnitTimeTrial_stmLaser,0,2);
rsUnitTimeTrial_stmLaserB_stmLaserZ=(rsUnitTimeTrial_stmLaser-repmat(meanPerUnit,[1 size(rsUnitTimeTrial_stmLaser,2)]))./repmat(stdPerUnit,[1 size(rsUnitTimeTrial_stmLaser,2)]);
unitTimeTrialBZ_stmLaser = reshape(rsUnitTimeTrial_stmLaserB_stmLaserZ, [size(unitTimeTrialB_stmLaser,1), size(unitTimeTrialB_stmLaser,2), size(unitTimeTrialB_stmLaser,3)]);

for t = 1:size(unitTimeTrialBZ_stmLaser,3) % trial-by-trial
    uttBZ_stmLaser = unitTimeTrialBZ_stmLaser(unitIdx,:,:); % unit indexing
    tmpNtj_stmLaser = projectMat'*uttBZ_stmLaser(:,:,t); % project the trial-averaged population activity of each fold to the projection matrix
    nTjCellAllTrials_stmLaser{1,1,t} = cell2mat(arrayfun(@(ROWIDX) conv(tmpNtj_stmLaser(ROWIDX,:),gaussianKernel2,'same'), (1:size(tmpNtj_stmLaser,1))', 'UniformOutput', false)); % smooth
    for uu = 1:size(uttBZ_stmLaser,1)
        tmpNtjEach_stmLaser = projectMat(uu,:)'*uttBZ_stmLaser(uu,:,t); % projection of each neuron
        nTjEachAllTrials_stmLaser{uu,1,t} = cell2mat(arrayfun(@(ROWIDX) conv(tmpNtjEach_stmLaser(ROWIDX,:),gaussianKernel2,'same'), (1:size(tmpNtjEach_stmLaser,1))', 'UniformOutput', false)); % smooth
    end
end
clearvars t

nTjC.stmLaser = nanmean(reshape([nTjCellAllTrials_stmLaser{:}],3,100,[]),3);

% trajectories of stmLaser trials per unit
for uu = 1:size(nTjEachAllTrials_stmLaser,1)
    prjKN_unit.stmLaser{uu,1} = cell2mat(nTjEachAllTrials_stmLaser(uu,:,:)); 
end

% trial-by-trial projection scores
for dd = 1:3
    prjKN_stmLaser{dd,1} = cell2mat(squeeze(cellfun(@(a) a(dd,:), nTjCellAllTrials_stmLaser, 'un', 0)));
end

save(fullfile(targDir.folder, targDir.name), 'prjKN_unit*', '-append')
end

