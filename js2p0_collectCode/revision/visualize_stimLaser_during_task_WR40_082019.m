% For visualization, this needs to be run first: js2p0_tbytSpkHandJsPreprocess_50ms_stim_parse(filePath)

filePath = '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles'; 
load(fullfile(filePath, strcat('js2p0_tbytSpkHandJsTrjBin_50ms_stimParse_', 'WR40_082019')),'ss', 'jkvt'); 
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'js2p0_behavior_stat_rez_inactivation.mat'), 'rezCol', 'rez', 'stat', 'trI') 
figSaveDir = '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles/Figure'; 

%% visualize representative trials
time_bins = (-975:50:3975)/1000;

% Create the meshgrid
[X, Y] = meshgrid(time_bins, 1:93);

t=75; 
figure; 
imagesc(X(1,:), Y(:,1), smooth2a(ss(t).utbCtxStimAlignZ, 0, 3)); 
colormap('hot'); 
clim([-1 3])

hold on;
% mark rStart
line([ss(t).rStartRelToStim./1000 ss(t).rStartRelToStim./1000], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
% mark laser on
line([0 0], ylim, 'Color', [135 206 235]./255, 'LineStyle', '--', 'LineWidth', 2);
% mark laser off
laserOff = (jkvt(t).stimLaserOff-jkvt(t).stimLaserOn)./1000; 
line([laserOff laserOff], ylim, 'Color', [135 206 235]./255, 'LineStyle', '--', 'LineWidth', 2);

% Label the axes
xlabel('Time (ms)');
set(gca, 'TickDir', 'out')
set(gca, 'XTick', -1:1:5)
%title('Smoothed Data Visualization');
hold off;

print(fullfile(figSaveDir, strcat('utbCtxStimAlignZ_', sprintf('trial#%d', t))), '-dpdf', '-bestfit', '-vector'); 





