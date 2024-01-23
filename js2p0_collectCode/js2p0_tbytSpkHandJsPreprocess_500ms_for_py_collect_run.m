filePath = {
    '/Volumes/dudmanlab/junchol/js2p0/WR37_022619/Matfiles', ... % Cg recording contra-Cg silencing (checked)
    '/Volumes/dudmanlab/junchol/js2p0/WR38_052119/Matfiles', ... % Dual recording without silencing (checked)
    '/Volumes/dudmanlab/junchol/js2p0/WR38_052219/Matfiles', ... % Dual recording without silencing (checked)
    '/Volumes/dudmanlab/junchol/js2p0/WR38_052319/Matfiles', ... % Cg recording contra-Cg silencing (checked)
    '/Volumes/dudmanlab/junchol/js2p0/WR38_052419/Matfiles', ... % Corticostriatal recording M1 silencing (checked)
    '/Volumes/dudmanlab/junchol/js2p0/WR39_100219/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/dudmanlab/junchol/js2p0/WR40_081919/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/dudmanlab/junchol/js2p0/WR40_082019/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/dudmanlab/junchol/js2p0/WR44_031020/Matfiles'};    % Dual recording with contra Cg delayed silencing (checked)

for f = 1:length(filePath) 
    js2p0_tbytSpkHandJsPreprocess_500ms_for_py(filePath{f})
end

dat_paths = [    '/Volumes/dudmanlab/junchol/js2p0/WR37_022619/Matfiles', 
    '/Volumes/dudmanlab/junchol/js2p0/WR38_052119/Matfiles', 
    '/Volumes/dudmanlab/junchol/js2p0/WR38_052219/Matfiles', 
    '/Volumes/dudmanlab/junchol/js2p0/WR38_052319/Matfiles', 
    '/Volumes/dudmanlab/junchol/js2p0/WR38_052419/Matfiles', 
    '/Volumes/dudmanlab/junchol/js2p0/WR39_100219/Matfiles', 
    '/Volumes/dudmanlab/junchol/js2p0/WR40_081919/Matfiles', 
    '/Volumes/dudmanlab/junchol/js2p0/WR40_082019/Matfiles', 
    '/Volumes/dudmanlab/junchol/js2p0/WR44_031020/Matfiles']