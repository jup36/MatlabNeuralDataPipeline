filePaths = {
    '/Volumes/Extreme SSD/js2p0/WR38_052119/Matfiles', ... % Dual recording without Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR38_052419/Matfiles', ... % Corticostriatal recording M1 silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR39_100219/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR44_031020/Matfiles'};    % Dual recording with contra Cg delayed silencing (checked)

mR2_tbyt = {}; 
mR2_tbytStim = {}; 

for f = 3:length(filePaths)
    filePath = GrabFiles_sort_trials('js2p0_tbytSpkHandJsTrjBin_50ms_', 0, filePaths(f));
    [~, rrrRezR2] = stimEffect_neuralTrj_cv(filePath{1}, 20, 10); % args: dimensions and folds 

    mR2_tbyt{f} = rrrRezR2.mR2_tbyt; 
    mR2_tbytStim{f} = rrrRezR2.mR2_tbytStim; 
end