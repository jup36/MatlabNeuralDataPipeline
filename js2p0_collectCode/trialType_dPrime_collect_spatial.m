
filePath = {'/Volumes/8TB/Junchol_Data/JS2p0/WR37_022119/Matfiles'...,
    '/Volumes/8TB/Junchol_Data/JS2p0/WR38_052219/Matfiles'...,
    '/Volumes/8TB/Junchol_Data/JS2p0/WR38_052419/Matfiles'...,
    '/Volumes/8TB/Junchol_Data/JS2p0/WR39_100219/Matfiles'...,
    '/Volumes/8TB/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
    '/Volumes/8TB/Junchol_Data/JS2p0/WR40_082019/Matfiles'...,
    '/Volumes/8TB/Junchol_Data/JS2p0/WR44_031020/Matfiles'};

% neuronal encoding types: lelo (1) >> le (2) >> lehi (3) >> hi (4) >> rihi (5) >> ri (6) >> rilo (7) >> lo (8) (counter-clockwise from left-low)
load(fullfile('/Volumes/8TB/Junchol_Data/JS2p0/collectData','dPrime_CtxStr_collectRez'),'dPrmRez')
isStr = [dPrmRez(:).isStr]; 

%% Striatum
strDepth = [dPrmRez(isStr).depth]; 

[strDepthN,strDepthE,strDepthB] = histcounts(strDepth, 2000:100:max(strDepth)); 
strEncodeType = cell2mat(cellfun(@(a) a(2), {dPrmRez(isStr).sigMaxCoord_2dProj}, 'un', 0)); % striatal neuronal encoding types (1-8 or NaN)

for dd = 1:max(strDepthB)
    depth_dPrmType_Str{dd} = histcounts(strEncodeType(strDepthB==dd),1:9)'; % count encoding types within each depth bin
end

str_depth_sumCell = cell2mat(cellfun(@(a) nansum(a(:)), depth_dPrmType_Str, 'un', 0)); % all striatal cells by depth
str_depth_sumType = nansum(cell2mat(depth_dPrmType_Str),2); 

% pure direction encoding neurons
str_dir_depth_count = cell2mat(cellfun(@(a) nansum(a([2,6],1)), depth_dPrmType_Str, 'un', 0)); % pure-direction encoding cells 
str_dir_depth_normBin = str_dir_depth_count./str_depth_sumCell; 
str_dir_depth_normType = str_dir_depth_count./nansum(str_dir_depth_count); 

% pure direction encoding neurons
str_load_depth_count = cell2mat(cellfun(@(a) sum(a([4,8],1)), depth_dPrmType_Str, 'un', 0)); % pure-load encoding cells
str_load_depth_normBin = str_load_depth_count./str_depth_sumCell; 
str_load_depth_normType = str_load_depth_count./nansum(str_load_depth_count); 

% mixed encoding neurons
str_mix_depth_count = cell2mat(cellfun(@(a) sum(a([1,3,5,7],1)), depth_dPrmType_Str, 'un', 0)); % mixed encoding cells
str_mix_depth_normBin = str_mix_depth_count./str_depth_sumCell; 
str_mix_depth_normType = str_mix_depth_count./nansum(str_mix_depth_count); 

% plot encoding type proportion per each depth bin
strX = strDepthE(1:end-1); 
figure; hold on; 
plot(strX,smooth2a(str_dir_depth_normBin,0,2))
plot(strX,smooth2a(str_load_depth_normBin,0,2))
plot(strX,smooth2a(str_mix_depth_normBin,0,2))
axis tight

% plot encoding type proportion per each depth bin normalized such that each encoding type's max to be 1
figure; hold on; 
plot(strX,smooth2a(str_dir_depth_normBin,0,2)./max(smooth2a(str_dir_depth_normBin,0,2)))
plot(strX,smooth2a(str_load_depth_normBin,0,2)./max(smooth2a(str_load_depth_normBin,0,2)))
plot(strX,smooth2a(str_mix_depth_normBin,0,2)./max(smooth2a(str_mix_depth_normBin,0,2)))
axis tight
set(gca,'tickDir','out')
print(fullfile('/Volumes/8TB/Junchol_Data/JS2p0/collectData/collectFigure','dPrm_encodingType_normProportion_Str'),'-dpdf','-painters')

% plot encoding type proportion per each encoding type 
figure; hold on; 
plot(strX,smooth2a(str_dir_depth_normType,0,2))
plot(strX,smooth2a(str_load_depth_normType,0,2))
plot(strX,smooth2a(str_mix_depth_normType,0,2))

%% Cortex
ctxDepth = [dPrmRez(~isStr).depth]; 

[ctxDepthN,ctxDepthE,ctxDepthB] = histcounts(ctxDepth, floor(min(ctxDepth)):100:1900); 
ctxEncodeType = cell2mat(cellfun(@(a) a(2), {dPrmRez(~isStr).sigMaxCoord_2dProj}, 'un', 0)); % cortical neuronal encoding types (1-8 or NaN)

for dd = 1:max(ctxDepthB)
    depth_dPrmType_Ctx{dd} = histcounts(ctxEncodeType(ctxDepthB==dd),1:9)'; % count encoding types within each depth bin
end

ctx_depth_sumCell = cell2mat(cellfun(@(a) nansum(a(:)), depth_dPrmType_Ctx, 'un', 0)); % all ctxiatal cells by depth
ctx_depth_sumType = nansum(cell2mat(depth_dPrmType_Ctx),2); 

% pure direction encoding neurons
ctx_dir_depth_count = cell2mat(cellfun(@(a) nansum(a([2,6],1)), depth_dPrmType_Ctx, 'un', 0)); % pure-direction encoding cells 
ctx_dir_depth_normBin = ctx_dir_depth_count./ctx_depth_sumCell; 
ctx_dir_depth_normType = ctx_dir_depth_count./nansum(ctx_dir_depth_count); 

% pure direction encoding neurons
ctx_load_depth_count = cell2mat(cellfun(@(a) sum(a([4,8],1)), depth_dPrmType_Ctx, 'un', 0)); % pure-load encoding cells
ctx_load_depth_normBin = ctx_load_depth_count./ctx_depth_sumCell; 
ctx_load_depth_normType = ctx_load_depth_count./nansum(ctx_load_depth_count); 

% mixed encoding neurons
ctx_mix_depth_count = cell2mat(cellfun(@(a) sum(a([1,3,5,7],1)), depth_dPrmType_Ctx, 'un', 0)); % mixed encoding cells
ctx_mix_depth_normBin = ctx_mix_depth_count./ctx_depth_sumCell; 
ctx_mix_depth_normType = ctx_mix_depth_count./nansum(ctx_mix_depth_count); 

% plot encoding type proportion per each depth bin
ctxX = ctxDepthE(1:end-1); 
figure; hold on; 
plot(ctxX,smooth2a(ctx_dir_depth_normBin,0,2))
plot(ctxX,smooth2a(ctx_load_depth_normBin,0,2))
plot(ctxX,smooth2a(ctx_mix_depth_normBin,0,2))

% plot encoding type proportion per each depth bin normalized such that each encoding type's max to be 1
figure; hold on; 
plot(ctxX,smooth2a(ctx_dir_depth_normBin,0,2)./max(smooth2a(ctx_dir_depth_normBin,0,2)))
plot(ctxX,smooth2a(ctx_load_depth_normBin,0,2)./max(smooth2a(ctx_load_depth_normBin,0,2)))
plot(ctxX,smooth2a(ctx_mix_depth_normBin,0,2)./max(smooth2a(ctx_mix_depth_normBin,0,2)))
axis tight
set(gca,'tickDir','out')
print(fullfile('/Volumes/8TB/Junchol_Data/JS2p0/collectData/collectFigure','dPrm_encodingType_normProportion_Ctx'),'-dpdf','-painters')

% plot encoding type proportion per each encoding type 
figure; hold on; 
plot(ctxX,smooth2a(ctx_dir_depth_normType,0,2))
plot(ctxX,smooth2a(ctx_load_depth_normType,0,2))
plot(ctxX,smooth2a(ctx_mix_depth_normType,0,2))


