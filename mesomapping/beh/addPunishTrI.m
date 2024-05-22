function tbytDat = addPunishTrI(tbytDat)
punishTrI = zeros(1, length(tbytDat));
if isfield(tbytDat, 'pos_rwd_tr')
    punishTrI = punishTrI + [tbytDat.neg_rwd_tr];
end

if isfield(tbytDat, 'pos_oms_tr')
    punishTrI = punishTrI + [tbytDat.neg_oms_tr];
end

if isfield(tbytDat, 'pos_pns_tr')
    punishTrI = punishTrI + [tbytDat.neg_pns_tr];
end

for i = 1:length(tbytDat)
    tbytDat(i).punishTrI = punishTrI(i);
end
end