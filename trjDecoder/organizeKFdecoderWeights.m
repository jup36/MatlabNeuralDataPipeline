% filePath = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR37_022119/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052219/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052419/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR39_100219/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'...,
%     '/Volumes/Beefcake/Junchol_Data/JS2p0/WR44_031020/Matfiles'};

function [ ctxRchW_posM, strRchW_posM, ctxRchW_velM, strRchW_velM, ctxPullW_posM, strPullW_posM, ctxPullW_velM, strPullW_velM ] = organizeKFdecoderWeights(filePath)

for f = 1:length(filePath)
    %% load files
    rDir = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPosVel_reach_*'));
    rS = load(fullfile(rDir.folder, rDir.name));
    pDir = dir(fullfile(filePath{f},'rezKFdecodeHTrjCtxStrPosVel_pull_*'));
    pS = load(fullfile(pDir.folder, pDir.name));
    
    bscDir = dir(fullfile(filePath{f},'binSpkCountSTRCTXWR*'));
    stc = load(fullfile(bscDir(1).folder,bscDir(1).name),'spkTimesCell');
    depths = cell2mat(cellfun(@(a) a(2), stc.spkTimesCell(4,:),'un', 0));
    
    % reach
    ctxRchW_posM{f,1} = collectAllWeights(rS.s.dat.params{1,7}.C_ctxVal,depths);
    strRchW_posM{f,1} = collectAllWeights(rS.s.dat.params{1,7}.C_strVal,depths);
    
    ctxRchW_velM{f,1} = collectAllWeights(rS.s.dat.params{1,8}.C_ctxVal,depths);
    strRchW_velM{f,1} = collectAllWeights(rS.s.dat.params{1,8}.C_strVal,depths);
    
    % pull
    ctxPullW_posM{f,1} = collectAllWeights(pS.s.dat.params{1,7}.C_ctxVal,depths);
    strPullW_posM{f,1} = collectAllWeights(pS.s.dat.params{1,7}.C_strVal,depths);
    
    ctxPullW_velM{f,1} = collectAllWeights(pS.s.dat.params{1,8}.C_ctxVal,depths);
    strPullW_velM{f,1} = collectAllWeights(pS.s.dat.params{1,8}.C_strVal,depths);
    
    % organize corr
    fprintf('processed file# %d\n', f)
end

    function [ meanW ] = collectAllWeights(weightMat,depths)
        nW = numel(weightMat);
        meanW = nanmean(cell2mat(reshape(weightMat,[1,1,nW])),3); % average across all trial resamples
        meanW(:,4) = depths'./1000; % add the depth in the last column to sort
        meanW(:,5) = 1:size(meanW,1);
        meanW = sortrows(meanW(~isnan(sum(meanW,2)),:),4); % sort weights by depths
    end

end
