filePaths = {
    '/Volumes/Extreme SSD/js2p0/WR38_052119/Matfiles', ... % Dual recording without Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR38_052419/Matfiles', ... % Corticostriatal recording M1 silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR39_100219/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR44_031020/Matfiles'};    % Dual recording with contra Cg delayed silencing (checked)

mR2.stim = {};
mR2.pStim = {};
mR2.noStim = {}; 

mR2.stim_tbyt = {};
mR2.pStim_tbyt = {};
mR2.noStim_tbyt = {}; 

% preprocessing
%for f = 1:length(filePaths)
%    js2p0_tbytSpkHandJsPreprocess_50ms_stimPstim(filePaths{f})
%end

% train RRR with movement-potent corticostriatal activity and test them
% for stim- and pseudo-stim aligned activity
for f = 1:length(filePaths)
    filePath = GrabFiles_sort_trials('js2p0_tbytSpkHandJsTrjBin_50ms_stimPstim_', 0, filePaths(f));
    [~, rrrRezR2] = stimEffect_neuralTrj_cv_stim_pstim(filePath{1}, 50, 10); % args: dimensions and folds

    % organize r2 results
    mR2.noStim_tbyt{f} = rrrRezR2.mR2_tbyt;
    mR2.stim_tbyt{f}   = rrrRezR2.mR2_tbytStim;
    mR2.pStim_tbyt{f}  = rrrRezR2.mR2_tbytPstim;
    
    mR2.noStim{f} = rrrRezR2.mR2; 
    mR2.stim{f} = rrrRezR2.mR2_stim; 
    mR2.pStim{f} = rrrRezR2.mR2_pStim;

    medR2.noStim_tbyt{f} = rrrRezR2.medR2_tbyt;
    medR2.stim_tbyt{f}   = rrrRezR2.medR2_tbytStim;
    medR2.pStim_tbyt{f}  = rrrRezR2.medR2_tbytPstim;
    
    medR2.noStim{f} = rrrRezR2.medR2; 
    medR2.stim{f} = rrrRezR2.medR2_stim; 
    medR2.pStim{f} = rrrRezR2.medR2_pStim;
end








