

filePath = {'D:\Junchol_Data\JS2p0\WR37_022119', ... % B Only
    'D:\Junchol_Data\JS2p0\WR37_022219', ... % B Only
    'D:\Junchol_Data\JS2p0\WR37_022619', ... % Cg recording contra-Cg silencing
    'D:\Junchol_Data\JS2p0\WR37_022719', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR38_052119', ... % Dual recording without silencing
    'D:\Junchol_Data\JS2p0\WR38_052219', ... % Dual recording without silencing
    'D:\Junchol_Data\JS2p0\WR38_052319', ... % Cg recording contra-Cg silencing
    'D:\Junchol_Data\JS2p0\WR38_052419', ... % Corticostriatal recording M1 silencing
    'D:\Junchol_Data\JS2p0\WR39_091019', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR39_091119', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR39_100219', ... % Dual recording with contra Cg silencing
    'D:\Junchol_Data\JS2p0\WR39_100319', ... % B Only contra-Cg delayed silencing
    'D:\Junchol_Data\JS2p0\WR40_081919', ... % Dual recording with contra Cg silencing
    'D:\Junchol_Data\JS2p0\WR40_082019', ... % Dual recording with contra Cg silencing
    'D:\Junchol_Data\JS2p0\WR44_031020'};    % Dual recording with contra Cg delayed silencing


for f = 1:length(filePath)
    run_reachPull_trialType_dPrime_vector_norm(filePath{f})
    fprintf('processed file #%d\n', f)
end