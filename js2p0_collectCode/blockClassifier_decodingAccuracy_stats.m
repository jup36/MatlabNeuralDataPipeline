% decoding accuracy reampled
MOp = [0.6111, 0.5357, 0.5172, 0.5714, 0.6538, 0.6190, 0.5357].*100;
dSTR = [0.5789, 0.6071, 0.4828, 0.5172, 0.6667, 0.5714, 0.5714].*100; 
MOs = [0.4516, 0.6667, 0.4643, 0.4062, 0.4074, 0.4762, 0.3929].*100; 
    
data = [MOp, dSTR, MOs]';
group = [repmat({'MOp'}, length(MOp), 1); 
         repmat({'dSTR'}, length(dSTR), 1); 
         repmat({'MOs'}, length(MOs), 1)];

[p, tbl, stats] = anova1(data, group, 'off');  % 'off' means no default ANOVA plot

figure;
[c, m, h, nms] = multcompare(stats, 'CType', 'bonferroni');
