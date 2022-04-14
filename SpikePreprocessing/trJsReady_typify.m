function [trJsReady] = trJsReady_typify(jkvt)
trJsReady_trI = find(~cellfun(@isempty, {jkvt.trJsReady}));
rStartToPull_trI = find(~cellfun(@isempty, {jkvt.rStartToPull}));
pullStarts_trI = find(~cellfun(@isempty, {jkvt.pullStarts}));
reward_trI = find(cellfun(@(a) a == 1, {jkvt.rewarded}));

trJsReady.rStartToPullIdx = ismember(trJsReady_trI, rStartToPull_trI);
trJsReady.pullStartsIdx = ismember(trJsReady_trI, pullStarts_trI);
trJsReady.rewardIdx = ismember(trJsReady_trI, reward_trI);
end
