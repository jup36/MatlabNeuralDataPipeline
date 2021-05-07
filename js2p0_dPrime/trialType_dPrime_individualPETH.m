filePath = {'/Volumes/8TB/Junchol_Data/JS2p0/WR37_022119/Matfiles'...,
            '/Volumes/8TB/Junchol_Data/JS2p0/WR38_052219/Matfiles'...,
            '/Volumes/8TB/Junchol_Data/JS2p0/WR38_052419/Matfiles'...,
            '/Volumes/8TB/Junchol_Data/JS2p0/WR39_100219/Matfiles'...,
            '/Volumes/8TB/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
            '/Volumes/8TB/Junchol_Data/JS2p0/WR40_082019/Matfiles'...,
            '/Volumes/8TB/Junchol_Data/JS2p0/WR44_031020/Matfiles'};
saveName = {'WR37_022119','WR38_052219','WR38_052419','WR39_100219','WR40_081919','WR40_082019','WR44_031020'};

load('/Volumes/8TB/Junchol_Data/JS2p0/collectData/dPrime_CtxStr_collectRez.mat')

for tt = 2:length(dPrmTtC)
    for cc = 1:length(dPrmTtC{1,tt})
        cellId_dPrmTtC = dPrmTtC{1,tt}(cc,3);
        
        cellOrgLabel = dPrmRez(cellId_dPrmTtC).cellId;
        
        if dPrmRez(cellId_dPrmTtC).isStr
            regionId = 'Str';
        else
            regionId = 'Ctx';
        end
        cellFileId = cellOrgLabel(1:11);
        cellId = str2double(cellOrgLabel(strfind(cellOrgLabel,'#')+1:end));
        fileId = cell2mat(cellfun(@(a) contains(a, cellFileId), saveName, 'un', 0));
        typeId = dPrmRez(cellId_dPrmTtC).sigMaxCoord_2dProj(end);
        
        spkDir = dir(fullfile(filePath{fileId},'binSpkCountSTRCTX*'));
        load(fullfile(spkDir(1).folder, spkDir(1).name),'rStartToPull','jkvt','spkTimesCell');
        S = rStartToPull;
        
        %% individual unit psth aligned to a task event
        % sort trials
        clearvars tq tqT ps psT tqpsT
        tq = round([jkvt.pull_torque]./10); % torque
        tqT = tq(S.trI)'; tqT(:,2) = 1:length(S.trI); tqT(:,3) = S.trI;  sortByTq = sortrows(tqT,1);
        ps = [jkvt.reachP1]; % position
        psT = ps(S.trI)'; psT(:,2) = 1:length(S.trI); psT(:,3) = S.trI;  sortByPs = sortrows(psT,1);
        tqpsT(:,1) = tqT(:,1)+psT(:,1); tqpsT(:,2) = 1:length(S.trI); tqpsT(:,3) = S.trI; sortByTqPs = sortrows(tqpsT,1); % torque position combo
        
        [cellI] = getUnitIdSfromSpkTimesCell(spkTimesCell, S, cellId);
        thisUnitSpkTimes = S.SpkTimes{cellI};
        %individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByTq, cellI, [3e3 2e3], [1e3 1e3]);
        %individualUnitPlotSortByType(filePath, thisUnitSpkTimes, sortByPs, cellI, [3e3 2e3], [1e3 1e3]);
        saveName1 = strcat(num2str(cellId_dPrmTtC),'_',sprintf('Type#%d',typeId),'_',sprintf('stcID#%d_sId#%d',cellId,cellI),'_',saveName{fileId},'_',regionId);
        individualUnitPlotSortByType_saveName(fullfile('/Volumes/8TB/Junchol_Data/JS2p0/collectData/collectFigure/IndividualNeuronPSTH',sprintf('Type%d',typeId)),...
            thisUnitSpkTimes, sortByTqPs, [3e3 2e3], [1e3 1e3], saveName1);
        close all;
    end
end