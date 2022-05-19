
filePaths = {'D:\Junchol_Data\JS2p0\WR37_022119', ... % B Only
    'D:\Junchol_Data\JS2p0\WR37_022219', ... % B Only
    'D:\Junchol_Data\JS2p0\WR37_022619', ... % Cg recording contra-Cg silencing
    'D:\Junchol_Data\JS2p0\WR37_022719', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR38_052119', ... % Dual recording without silencing
    'D:\Junchol_Data\JS2p0\WR38_052219', ... % Dual recording without silencing
    'D:\Junchol_Data\JS2p0\WR38_052319', ... % Cg recording contra-Cg silencing
    'D:\Junchol_Data\JS2p0\WR38_052419', ... % Corticostriatal recording M1 silencing
    'D:\Junchol_Data\JS2p0\WR39_091019', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR39_091119', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR39_100219', ... % Dual recording with contra Cg silencing
    'D:\Junchol_Data\JS2p0\WR39_100319', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR40_081919', ... % Dual recording with contra Cg silencing
    'D:\Junchol_Data\JS2p0\WR40_082019', ... % Dual recording with contra Cg silencing
    'D:\Junchol_Data\JS2p0\WR44_031020'};    % Dual recording with contra Cg delayed silencing

load(fullfile('D:\Junchol_Data\JS2p0\collectData','dPrime_CtxStr_vec_norm_collectRez'),'dPrmRez')

%% choose what to draw per block
% CTX
for bb = 1:8
    regionI = cell2mat(cellfun(@(a) a==0, {dPrmRez.isStr}, 'un', 0));
    blockI = cell2mat(cellfun(@(a) a(2)==bb, {dPrmRez.maxCoord_2dProj}, 'un', 0));
    dPrm = dPrmRez(blockI & regionI);
    
    vec_norm = cell2mat(cellfun(@(a) a(1), {dPrm.maxCoord_2dProj}, 'un', 0))';
    vec_norm(:,2) = 1:length(vec_norm);
    
    sorted = sortrows(vec_norm, -1);
    
    dPrm_plot = dPrm(sorted(1:10,2)); % this is to be plotted
    
    rank = 1; 
    for cc = 1:length(dPrm_plot)
        filePath = filePaths{cell2mat(cellfun(@(a) contains(a, dPrm_plot(cc).cellId(1:11)), filePaths, 'un', 0))};
        cd(filePath)
        fileDir_ss = dir('**/*js2p0_tbytSpkHandJsTrjBin_WR*.mat');
        load(fullfile(fileDir_ss(1).folder, fileDir_ss(1).name), 'ss');
        fileDir_bsc = dir('**/*binSpkCountSTRCTXWR*.mat');
        load(fullfile(fileDir_bsc(1).folder, fileDir_bsc(1).name), 'rStartToPull', 'jkvt', 'spkTimesCell')
        
        S = rStartToPull;
        cellI = str2double(dPrm_plot(cc).cellId(strfind(dPrm_plot(cc).cellId,'#')+1:end)); 
        [cellId] = getUnitIdSfromSpkTimesCell(spkTimesCell, S, cellI);
        thisUnitSpkTimes = S.SpkTimes{cellId};
        % sort trials
        b_id = [ss.blNumber]'; 
        b_id = b_id(S.trI); 
        
        b_type = zeros(length(b_id),1); 
        for t = 1:length(b_id)
            if b_id(t) == 1 || b_id(t) == 5 
               b_type(t)= 1; 
            elseif b_id(t) == 2 || b_id(t) == 6 
               b_type(t)= 2;  
            elseif b_id(t) == 3 || b_id(t) == 7 
               b_type(t)= 3;   
            elseif b_id(t) == 4 || b_id(t) == 8 
               b_type(t)= 4;   
            end 
        end
        b_type(:,2) = 1:length(b_id); 
        % plot and save
        individualUnitPlotSortByType_figSave_each(filePath, dPrm_plot(cc), thisUnitSpkTimes, b_type, cellI, [3e3 2e3], [1e3 1e3], rank);
        clearvars tqpsT sortByTqPs
        fprintf('processed rank # %d\n', rank) % report unit progression
        rank = rank + 1; 
    end
    fprintf('processed block # %d\n', bb) % report unit progression
end

% STR
for bb = 1:8
    regionI = cell2mat(cellfun(@(a) a==1, {dPrmRez.isStr}, 'un', 0));
    blockI = cell2mat(cellfun(@(a) a(2)==bb, {dPrmRez.maxCoord_2dProj}, 'un', 0));
    dPrm = dPrmRez(blockI & regionI);
    
    vec_norm = cell2mat(cellfun(@(a) a(1), {dPrm.maxCoord_2dProj}, 'un', 0))';
    vec_norm(:,2) = 1:length(vec_norm);
    
    sorted = sortrows(vec_norm, -1);
    
    dPrm_plot = dPrm(sorted(1:10,2)); % this is to be plotted
    
    rank = 1; 
    for cc = 1:length(dPrm_plot)
        filePath = filePaths{cell2mat(cellfun(@(a) contains(a, dPrm_plot(cc).cellId(1:11)), filePaths, 'un', 0))};
        cd(filePath)
        fileDir_ss = dir('**/*js2p0_tbytSpkHandJsTrjBin_WR*.mat');
        load(fullfile(fileDir_ss(1).folder, fileDir_ss(1).name), 'ss');
        fileDir_bsc = dir('**/*binSpkCountSTRCTXWR*.mat');
        load(fullfile(fileDir_bsc(1).folder, fileDir_bsc(1).name), 'rStartToPull', 'jkvt', 'spkTimesCell')
        
        S = rStartToPull;
        cellI = str2double(dPrm_plot(cc).cellId(strfind(dPrm_plot(cc).cellId,'#')+1:end)); 
        [cellId] = getUnitIdSfromSpkTimesCell(spkTimesCell, S, cellI);
        thisUnitSpkTimes = S.SpkTimes{cellId};
        % sort trials
        b_id = [ss.blNumber]'; 
        b_id = b_id(S.trI); 
        
        b_type = zeros(length(b_id),1); 
        for t = 1:length(b_id)
            if b_id(t) == 1 || b_id(t) == 5 
               b_type(t)= 1; 
            elseif b_id(t) == 2 || b_id(t) == 6 
               b_type(t)= 2;  
            elseif b_id(t) == 3 || b_id(t) == 7 
               b_type(t)= 3;   
            elseif b_id(t) == 4 || b_id(t) == 8 
               b_type(t)= 4;   
            end 
        end
        b_type(:,2) = 1:length(b_id); 
        % plot and save
        individualUnitPlotSortByType_figSave_each(filePath, dPrm_plot(cc), thisUnitSpkTimes, b_type, cellI, [3e3 2e3], [1e3 1e3], rank);
        clearvars tqpsT sortByTqPs
        fprintf('processed rank # %d\n', rank) % report unit progression
        rank = rank + 1; 
    end
    fprintf('processed block # %d\n', bb) % report unit progression
end



%%
function [unitIdS] = getUnitIdSfromSpkTimesCell(spkTimesCellIn, sIn, unitIdC)
%unit Ids might not match between 'spkTimesCell' and the structure 'S' (e.g. S.unitTimeTrial)
% this function identifies the unit from 'S' in 'spkTimesCell' and outputs
% it as 'unitidSTC'

geomI = cell2mat(cellfun(@(a) isequal(a,spkTimesCellIn{4,unitIdC}), sIn.geometry, 'un', 0));
wfI = cell2mat(cellfun(@(a) isequal(a,spkTimesCellIn{6,unitIdC}), sIn.meanWF, 'un', 0));

unitIdS = find(geomI & wfI);

end


function individualUnitPlotSortByType_figSave_each(filePath, dPrm, dataS, type, unit, psthWin, manualX, rank_within_block)
%This function is to generate an individual unit plot whose trials are
% sorted based on a behavioral variable by the # of fold.

% put data into a cell - each fold as one cell element

figSavePath = 'D:\Junchol_Data\JS2p0\collectData\collectFigure\IndividualNeuronPSTH';
wrI = strfind(filePath, 'WR');
fileName = filePath(wrI:wrI+10);

block = num2str(dPrm.maxCoord_2dProj(2));

types = unique(type(:,1)); % trial p.Results.trialFolds
typeCell = cell(1,length(types));

% get the trial-averaged neural population trajectories of each fold rank-ordered by a movement variable
for f = 1:length(types)
    typeCell{1,f} = dataS(type(type(:,1)==types(f),2)); % take the trials of the current fold sorted by beh
end
clearvars f

unitIdFmt = 'Unit%d';
rankFmt = 'Rank%d';
if dPrm.isStr==0
    saveName = strcat(fileName,'_Ctx','_block', block, '_',sprintf(unitIdFmt, unit), '_', sprintf(rankFmt, rank_within_block));
else
    saveName = strcat(fileName,'_Str','_block', block, '_',sprintf(unitIdFmt, unit), '_', sprintf(rankFmt, rank_within_block));
end

h = spikeRasterGrammSortedFolds( psthWin, manualX, typeCell ); % inputs; psthWin, manualX, typeCell
pbaspect([1 1 1])

print(h,fullfile(figSavePath, strcat('Type',block), saveName), '-dpdf','-bestfit')
close all; 
end
