 filePaths = {
    '/Volumes/Extreme SSD/js2p0/WR38_052119/Matfiles', ... % Dual recording without Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR38_052419/Matfiles', ... % Corticostriatal recording M1 silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR39_100219/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR44_031020/Matfiles'};    % Dual recording with contra Cg delayed silencing (checked)

for f = 1:length(filePaths) 
    js2p0_tbytSpkHandJsPreprocess_50ms_stimPstimPrepExt(filePaths{f})
    fprintf('Completed preprocessing file #%d\n', f)
end