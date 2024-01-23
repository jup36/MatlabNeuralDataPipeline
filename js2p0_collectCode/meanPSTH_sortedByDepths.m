
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


binUnitTimeCat = [];
depth = [];

for f = 1:length(filePath)
    cd(filePath{f})
    
    spkDir = dir('**/*binSpkCountSTRCTX*.mat');
    if ~isempty(spkDir)
        load(fullfile(spkDir(1).folder, spkDir(1).name), 'rStartToPull')
        
        unitTimeTrMean = nanmean(rStartToPull.unitTimeTrial, 3);
        binUnitTime = bin1msSpkCountMat(unitTimeTrMean, 50, 50);
        
        binUnitTimeCat = [binUnitTimeCat; binUnitTime];
        depth = [depth; max(0, cell2mat(cellfun(@(a) a(2), rStartToPull.geometry, 'un', 0))-100)];
    end
    fprintf('processed file # %d\n', f) % report unit progression
end

[~, depthEdge, depthBin] = histcounts(depth,0:50:4000); 


depthZ = nan(length(depthEdge)-1, size(binUnitTime,2)); 
for dd = 1:length(depthEdge)-1
    depthZ(dd,:) = nanmean(binUnitTimeCat(depthBin==dd,:));
end

cmap  = TNC_CreateRBColormapJP(200,'rb'); % generate a colormap for imagesc psth

x = -2925:50:2025; 
y = depthEdge(4:end-1)/1000; 
imagesc(x, y, smooth2a(depthZ(4:end,:),1,1))
colorbar % units with a significant inhibitory tagging effect
colormap(cmap);
caxis([-1, 1.5])
xlim([-1000 1000])
set(gca,'TickDir','out')
print(fullfile('D:\Junchol_Data\JS2p0\collectData\collectFigure', 'meanPSTH_sortedByDepth'), '-dpdf', '-painters')



