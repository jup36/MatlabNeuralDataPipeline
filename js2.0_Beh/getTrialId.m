function [combinedIds, valTrialId, valBlockId] = getTrialId(fileName, trialI, blockI, datCell, refDat)
% datCell = {ss.rchAngDeg};
% trialI = trI.spNoStimI{f};
% blockI = rez.b_id{f};
% refDat = rez.absRchAng{f};

valDatI = ~cellfun(@isempty, datCell); 

valTrialI = valDatI(:) & trialI(:);
valTrialId = find(valTrialI);
assert(sum(valTrialI)==length(refDat))
valBlockId = blockI(valTrialI);
combinedIds = cell(sum(valTrialI), 1);

for i = 1:length(valBlockId)
    switch valBlockId(i)
        case {1, 5}
            combinedIds{i} = strcat(fileName, '_', sprintf("tr%d_bl%d", valTrialId(i), valBlockId(i)), '_', 'lelo');
        case {2, 6}
            combinedIds{i} = strcat(fileName, '_', sprintf("tr%d_bl%d", valTrialId(i), valBlockId(i)), '_', 'lehi');
        case {3, 7}
            combinedIds{i} = strcat(fileName, '_', sprintf("tr%d_bl%d", valTrialId(i), valBlockId(i)), '_', 'rilo');
        case {4, 8}
            combinedIds{i} = strcat(fileName, '_', sprintf("tr%d_bl%d", valTrialId(i), valBlockId(i)), '_', 'rihi');
    end
end
end