function dat = tbytDatRewardPunishI(dat)
rewardTrI = zeros(1, length(dat));
punishTrI = zeros(1, length(dat));

if isfield(dat, 'pos_rwd_tr')
    rewardTrI = rewardTrI+[dat.pos_rwd_tr];
    if isfield(dat, 'pos_oms_tr')
        rewardTrI = rewardTrI+[dat.pos_oms_tr];
        if isfield(dat, 'pos_pns_tr')
            rewardTrI = rewardTrI+[dat.pos_pns_tr];
        end
    end
end

if isfield(dat, 'neg_rwd_tr')
    punishTrI = punishTrI+[dat.neg_rwd_tr];
    if isfield(dat, 'neg_oms_tr')
        punishTrI = punishTrI+[dat.neg_oms_tr];
        if isfield(dat, 'neg_pns_tr')
            punishTrI = punishTrI+[dat.neg_pns_tr];
        end
    end
end

for i = 1:length(dat)
    dat(i).rewardTrI = rewardTrI(i);
    dat(i).punishTrI = punishTrI(i);

end
end
