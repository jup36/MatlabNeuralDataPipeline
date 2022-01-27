%% Helper functions
function plot_handJsTrjXY1_clockwiseRotation(filePath, varargin)

%filePath = fullfile('D:\Junchol_Data\JS2p0\WR40_081919\Matfiles','js2p0_tbytSpkHandJsTrjBin_WR40_081919.mat'); 
load(fullfile(filePath), 'ss', 'jkvt')
valTrs = cell2mat(cellfun(@(a) ~isempty(a), {ss(:).hTrjB}, 'un', 0));
rwdTrI = [jkvt(:).rewarded]; 










end