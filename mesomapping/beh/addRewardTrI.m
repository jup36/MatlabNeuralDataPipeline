function tbytDat = addRewardTrI(tbytDat)
rewardTrI = zeros(1, length(tbytDat));
if isfield(tbytDat, 'pos_rwd_tr')
    rewardTrI = rewardTrI + [tbytDat.pos_rwd_tr];
end

if isfield(tbytDat, 'pos_oms_tr')
    rewardTrI = rewardTrI + [tbytDat.pos_oms_tr];
end

if isfield(tbytDat, 'pos_pns_tr')
    rewardTrI = rewardTrI + [tbytDat.pos_pns_tr];
end

for i = 1:length(tbytDat)
    tbytDat(i).rewardTrI = rewardTrI(i);
end
end