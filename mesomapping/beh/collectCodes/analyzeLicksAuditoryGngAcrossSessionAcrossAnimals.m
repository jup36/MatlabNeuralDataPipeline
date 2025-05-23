
filePaths = {'/Volumes/buschman/Rodent Data/dualImaging_parkj/m1237_GCAMP', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1092_jRGECO', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1094_jRGECO', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1044_jRGECO_GRABda', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1048_jRGECO_GRABda', ...
             '/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda', ...
             };

dPrmCollectC = cell(length(filePaths), 1); 
hitCollectC = cell(length(filePaths), 1); 
crCollectC = cell(length(filePaths), 1); 

dPrmTcollectC = cell(length(filePaths), 1);
hitTcollectC = cell(length(filePaths), 1); 
crTcollectC = cell(length(filePaths), 1); 


for f = 1:length(filePaths)
    currPath = GrabFiles_sort_trials('_blockWise_behavior', 0, filePaths(f)); 
    load(currPath{1}, 'rezB', 'rezT'); 

    dPrmCollectC{f, 1} = rezB.dPrmC(:, 2)'; 
    hitCollectC{f, 1} = rezB.hitC(:, 2)'; 
    crCollectC{f, 1} = rezB.crC(:, 2)'; 

    dPrmTcollectC{f, 1} = rezT.dPrmTot'; 
    hitTcollectC{f, 1} = rezT.hitRateTot'; 
    crTcollectC{f, 1} = rezT.crRateTot'; 

    clearvars rezB rezT
end

%% plot
cMat = slanCM('vivid', length(filePaths)); 

% d prime all mice blocks
h_dPrm_all_blocks = plotLearningCurveCell(dPrmCollectC, cMat); 
ylabel('d prime')
print(h_dPrm_all_blocks, fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure', 'dPrime_learning_curves_collect_allAnimalsBlocks'), ...
    '-dpdf', '-vector', '-bestfit')

% d prime all mice 
h_dPrm_all = plotLearningCurve(dPrmTcollectC, cMat); 
ylabel('d prime')
print(h_dPrm_all, fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure', ...
    'dPrm_learning_curves_collect_allAnimals'), '-dpdf', '-vector', '-bestfit')

% d prime learners
h_dPrm_learners = plotLearningCurve(dPrmTcollectC(1:5,:), cMat(1:5,:)); 
ylabel('d prime')
print(h_dPrm_learners, fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure', ...
    'dPrm_learning_curves_collect_learners'), '-dpdf', '-vector', '-bestfit')

% hit rate all mice
h_hit_all = plotLearningCurve(hitTcollectC, cMat); 
ylabel('Hit Rate')
print(h_hit_all, fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure', ...
    'hitRate_learning_curves_collect_allAnimals'), '-dpdf', '-vector', '-bestfit')

% hit rate learners
h_hit_learners = plotLearningCurve(hitTcollectC(1:5,:), cMat(1:5,:)); 
ylabel('Hit Rate')
ylim([0 1])
print(h_hit_learners, fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure', ...
    'hitRate_learning_curves_collect_learners'), '-dpdf', '-vector', '-bestfit')

% CR rate all mice
h_CR_all = plotLearningCurve(crTcollectC, cMat); 
ylabel('CR Rate')
print(h_CR_all, fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure', ...
    'crRate_learning_curves_collect_allAnimals'), '-dpdf', '-vector', '-bestfit')

% CR rate learners
h_CR_learners = plotLearningCurve(crTcollectC(1:5,:), cMat(1:5,:)); 
ylabel('CR Rate')
ylim([0 1])
print(h_CR_learners, fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure', ...
    'crRate_learning_curves_collect_learners'), '-dpdf', '-vector', '-bestfit')