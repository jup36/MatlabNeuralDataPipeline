    function [lat] = latencyOfFirstLicks(tbytDat, trRwdI, trPnsI, cueWin)

    % First lick latency across blocks
    blocks = divideTrials(length(tbytDat), 3, 10, 10);

    lickTimeC = cellfun(@(a, b) a-b, {tbytDat.Lick}, {tbytDat.evtOn}, 'un', 0);

    lickTimeCueC = cellfun(@(a) a(a>=cueWin(1) & a<=cueWin(2)), lickTimeC, 'UniformOutput', false);

    lickCounts = cell2mat(cellfun(@length, lickTimeCueC, 'UniformOutput', false));
    valTrI = lickCounts > 0;
    valRwdTrI = valTrI(:) & trRwdI(:);

    % First lick latency
    lat.allFstLat = cell(1, length(tbytDat)); lat.allFstLat(:) = {NaN};
    lat.rwdFstLat = cell(1, length(tbytDat)); lat.rwdFstLat(:) = {NaN};
    lat.allFstLat(valTrI) = cellfun(@(a) a(1), lickTimeCueC(valTrI), 'UniformOutput', false);
    lat.rwdFstLat(valRwdTrI) = cellfun(@(a) a(1), lickTimeCueC(valRwdTrI), 'UniformOutput', false);

    % First lick latency counting miss trials as cue duration 
    lat.allFstLatMiss = cell(1, length(tbytDat));
    lat.allFstLatMiss(valTrI) = lat.allFstLat(valTrI); 
    lat.allFstLatMiss(~valTrI) = {max(cueWin)};  

    lat.rwdFstLatMiss = cell(1, length(tbytDat)); lat.rwdFstLatMiss(:) = {NaN};
    lat.rwdFstLatMiss(valRwdTrI) = cellfun(@(a) a(1), lickTimeCueC(valRwdTrI), 'UniformOutput', false);
    lat.rwdFstLatMiss(~valRwdTrI & trRwdI) = {max(cueWin)}; 

    % across blocks
    allFstLatBlocks = cellfun(@(a) cell2mat(lat.allFstLat(a)), blocks, 'UniformOutput', false);
    lat.allFstLatBlockMean = cell2mat(cellfun(@nanmean, allFstLatBlocks, 'UniformOutput', false));

    rwdFstLatBlocks = cellfun(@(a) cell2mat(lat.rwdFstLat(a)), blocks, 'UniformOutput', false);
    lat.rwdFstLatBlockMean = cell2mat(cellfun(@nanmean, rwdFstLatBlocks, 'UniformOutput', false));
    
    allFstLatMissBlocks = cellfun(@(a) cell2mat(lat.allFstLatMiss(a)), blocks, 'UniformOutput', false);
    lat.allFstLatMissBlockMean = cell2mat(cellfun(@nanmean, allFstLatMissBlocks, 'UniformOutput', false));

    rwdFstLatMissBlocks = cellfun(@(a) cell2mat(lat.rwdFstLatMiss(a)), blocks, 'UniformOutput', false);
    lat.rwdFstLatMissBlockMean = cell2mat(cellfun(@nanmean, rwdFstLatMissBlocks, 'UniformOutput', false));

    % for No-Go trials
    if sum(trPnsI)>5 % if there are No-Go trials
        valPnsTrI = valTrI(:) & trPnsI(:);
        lat.pnsFstLat = cell(1, length(tbytDat)); lat.pnsFstLat(:) = {NaN};
        lat.pnsFstLat(valPnsTrI) = cellfun(@(a) a(1), lickTimeCueC(valPnsTrI), 'UniformOutput', false);
        
        lat.pnsFstLatMiss = cell(1, length(tbytDat)); lat.pnsFstLatMiss(:) = {NaN};
        lat.pnsFstLatMiss(valPnsTrI) = cellfun(@(a) a(1), lickTimeCueC(valPnsTrI), 'UniformOutput', false);
        lat.pnsFstLatMiss(~valPnsTrI & trPnsI) = {max(cueWin)}; 

        % across blocks
        pnsFstLatBlocks = cellfun(@(a) cell2mat(lat.pnsFstLat(a)), blocks, 'UniformOutput', false);
        lat.pnsFstLatBlockMean = cell2mat(cellfun(@nanmean, pnsFstLatBlocks, 'UniformOutput', false));

        pnsFstLatMissBlocks = cellfun(@(a) cell2mat(lat.pnsFstLatMiss(a)), blocks, 'UniformOutput', false);
        lat.pnsFstLatMissBlockMean = cell2mat(cellfun(@nanmean, pnsFstLatMissBlocks, 'UniformOutput', false));
    end

    end