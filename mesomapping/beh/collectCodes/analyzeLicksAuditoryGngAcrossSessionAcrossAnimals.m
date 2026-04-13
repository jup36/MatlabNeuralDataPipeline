
filePaths = {compatiblepath('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1237_GCAMP'), ...
             compatiblepath('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1092_jRGECO'), ...
             compatiblepath('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1094_jRGECO'), ...
             compatiblepath('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1044_jRGECO_GRABda'), ...
             compatiblepath('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1045_jRGECO_GRABda'), ...
             compatiblepath('Z:\Rodent Data\dualImaging_parkj\m1613_jRGECO_GRABda'), ...
             compatiblepath('Z:\Rodent Data\dualImaging_parkj\m1859_jRGECO_GRABda'), ...
             compatiblepath('Z:\Rodent Data\dualImaging_parkj\m1873_jRGECO_GRABda'), ...
             compatiblepath('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1048_jRGECO_GRABda'), ...
             compatiblepath('/Volumes/buschman/Rodent Data/dualImaging_parkj/m1049_jRGECO_GRABda'), ...
             };

figSaveDir = compatiblepath('/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure'); 

dPrmCollectC = cell(length(filePaths), 1); 
hitCollectC = cell(length(filePaths), 1); 
crCollectC = cell(length(filePaths), 1); 

mIdC = cell(length(filePaths), 1); 
dPrmTcollectC = cell(length(filePaths), 1);
hitTcollectC = cell(length(filePaths), 1); 
crTcollectC = cell(length(filePaths), 1); 

for f = 1:length(filePaths)

    % mouse ID
    mIdC{f, 1} = regexp(filePaths{f}, 'm\d{3,5}', 'match', 'once');

    % load blockwise + total behavior
    currPath = GrabFiles_sort_trials('_blockWise_behavior', 0, filePaths(f));
    load(currPath{1}, 'rezB', 'rezT');

    % ---------- blockwise ----------
    % d'
    dPrmCollectC{f, 1} = rezB.dPrmC(:, 1)';  % headers
    dPrmCollectC{f, 2} = rezB.dPrmC(:, 2)';  % d' values

    % hit rate
    hitCollectC{f, 1} = rezB.hitC(:, 1)';    % headers
    hitCollectC{f, 2} = rezB.hitC(:, 2)';    % hit rates

    % correct rejection rate
    crCollectC{f, 1} = rezB.crC(:, 1)';      % headers
    crCollectC{f, 2} = rezB.crC(:, 2)';      % CR rates

    % ---------- total (session-averaged) ----------
    % mean d'
    dPrmTcollectC{f, 1} = rezB.dPrmC(:, 1)'; % headers
    dPrmTcollectC{f, 2} = rezT.dPrmTot';     % mean d'

    % mean hit rate
    hitTcollectC{f, 1} = rezB.hitC(:, 1)';   % headers
    hitTcollectC{f, 2} = rezT.hitRateTot';   % mean hit rate

    % mean CR rate
    crTcollectC{f, 1} = rezB.crC(:, 1)';     % headers
    crTcollectC{f, 2} = rezT.crRateTot';     % mean CR rate

    clearvars rezB rezT
end

%% learning rate defined as # of sessions to achieve the criterion (d'=1.5)
dPrmC = cell(size(dPrmCollectC, 1), 1); 

for j = 1:size(dPrmCollectC, 1)
    % dPrm
    dPrmC{j, 1} = dPrmCollectC{j, 1}; 
    dPrmC{j, 2} = cellfun(@(a) mean(a, 'omitnan'), dPrmCollectC{j, 2}); 
    [~, dPrmC{j, 3}] = find(dPrmC{j, 2}>1.5, 1, 'first'); 
end

%% save
save(compatiblepath(fullfile('Z:\Rodent Data\dualImaging_parkj\collectData', 'collect_behavior_dPrm_hit_cr_rates')), ...
    'dPrmCollectC', 'dPrmTcollectC', 'hitCollectC', 'crTcollectC', ... 
    'hitTcollectC', 'crTcollectC', 'dPrmC'); 

%% plot
cMat = slanCM('vivid', length(filePaths)); 

% d prime all mice blocks
h_dPrm_all_blocks = plotLearningCurveCell(dPrmCollectC, cMat); 
ylabel('d prime')
print(h_dPrm_all_blocks, fullfile(figSaveDir, 'dPrime_learning_curves_collect_allAnimalsBlocks'), ...
    '-dpdf', '-vector', '-bestfit')

% d prime all mice 
h_dPrm_all = plotLearningCurve(dPrmTcollectC, cMat); 
ylabel('d prime')
print(h_dPrm_all, fullfile(figSaveDir, ...
    'dPrm_learning_curves_collect_allAnimals'), '-dpdf', '-vector', '-bestfit')

% d prime each mouse
plotLearningCurveEachMouse(dPrmTcollectC, cMat, figSaveDir, 'learningCurveDprm'); 

% d prime learners
h_dPrm_learners = plotLearningCurve(dPrmTcollectC(1:8,:), cMat(1:8,:)); 
ylabel('d prime')
print(h_dPrm_learners, fullfile(figSaveDir, ...
    'dPrm_learning_curves_collect_learners'), '-dpdf', '-vector', '-bestfit')

% hit rate all mice
h_hit_all = plotLearningCurve(hitTcollectC, cMat); 
ylabel('Hit Rate')
print(h_hit_all, fullfile(figSaveDir, ...
    'hitRate_learning_curves_collect_allAnimals'), '-dpdf', '-vector', '-bestfit')

% hit rate learners
h_hit_learners = plotLearningCurve(hitTcollectC(1:8,:), cMat(1:8,:)); 
ylabel('Hit Rate')
ylim([0 1])
print(h_hit_learners, fullfile(figSaveDir, ...
    'hitRate_learning_curves_collect_learners'), '-dpdf', '-vector', '-bestfit')

% CR rate all mice
h_CR_all = plotLearningCurve(crTcollectC, cMat); 
ylabel('CR Rate')
print(h_CR_all, fullfile(figSaveDir, ...
    'crRate_learning_curves_collect_allAnimals'), '-dpdf', '-vector', '-bestfit')

% CR rate learners
h_CR_learners = plotLearningCurve(crTcollectC(1:8,:), cMat(1:8,:)); 
ylabel('CR Rate')
ylim([0 1])
print(h_CR_learners, fullfile(figSaveDir, ...
    'crRate_learning_curves_collect_learners'), '-dpdf', '-vector', '-bestfit')