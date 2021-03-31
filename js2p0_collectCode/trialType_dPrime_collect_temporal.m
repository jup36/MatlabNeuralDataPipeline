
filePath = {'/Volumes/Beefcake/Junchol_Data/JS2p0/WR37_022119/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052219/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR38_052419/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR39_100219/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_082019/Matfiles'...,
    '/Volumes/Beefcake/Junchol_Data/JS2p0/WR44_031020/Matfiles'};

% completed list:
load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','a2dColorMap.mat'),'colormap2D') % 2d colorMap for the scatter plot
numColors = size(colormap2D,1);
rangeColor = linspace(-2,2,numColors);

% lelo>> le >> lehi >> hi >> rihi >> ri >> rilo >> lo (counter-clockwise from left-low)
q2=sqrt(2);
projMatX = [-1 -q2 -1 0  1 q2  1   0];
projMatY = [-1  0   1 q2 1  0 -1 -q2];

projMat  = [projMatX; projMatY];
timeX = -1000:20:1000; 
timeXI = timeX>=0;% & timeX<500; 

load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','dPrime_CtxStr_collectRez'),'dPrmRez')

%% get dPrm temporal trajectories
for tt = 1:length(projMatX)
    [~,colorXI] = min(abs(projMat(1,tt)-rangeColor));
    [~,colorYI] = min(abs(projMat(2,tt)-rangeColor));
    colorC{tt} = [colormap2D(colorXI, colorYI, 1), colormap2D(colorXI, colorYI, 2), colormap2D(colorXI, colorYI, 3)];
end

for cc = 1:length(dPrmRez)
    tType = dPrmRez(cc).sigMaxCoord_2dProj(2); 
    
    if ~isempty(dPrmRez(cc).s) && ~isnan(tType)
       dPrmTrjC{cc,tType} = (dPrmRez(cc).s.trj2d*projMat(:,tType))'; 
       dPrmTrjC_distZ{cc,tType} = dPrmRez(cc).s.trjDistZero'; 
    end
end 

% figure; hold on; 
% for tt = 1:8
%     [mdPrm,~,sdPrm] = meanstdsem(cell2mat(dPrmTrjC(:,tt))); 
%     plot(mdPrm./max(mdPrm), 'Color', colorC{tt},'LineWidth',2)
% end

%% 
for tt = 1:8
    dPrmTimeTrjThisType = cell2mat(dPrmTrjC([dPrmRez(:).isStr]',tt)); % striatum 
    for cc = 1:size(dPrmTimeTrjThisType,1)
        dPrmTrjC_norm{tt}(cc,:) = dPrmTimeTrjThisType(cc,:)./max(dPrmTimeTrjThisType(cc,:),[],2); 
    end 
end
% plot pure load
figure; hold on; 
%plot(mean(mdPrm_str([4,8],:))./max(mean(mdPrm_str([4,8],:))))
%plot(mean(mdPrm_str([4,8],:)))
plot(mdPrm_str(8,:)./max(mdPrm_str(8,:)))
plot(mdPrm_str(4,:)./max(mdPrm_str(4,:)))

% plot pure direction 
%plot(mean(mdPrm([2,6],:))./max(mean(mdPrm([2,6],:))))
plot(mean(mdPrm_str([2,6],:)))

% plot pure mixed 
%plot(mean(mdPrm([1,3,5,7],:))./max(mean(mdPrm([1,3,5,7],:))))
plot(mean(mdPrm_str([1,3,5,7],:)))
axis tight
set(gca,'tickDir','out')
print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData/collectFigure','dPrm_neuronType_Str'),'-dpdf','-painters')

%%
for tt = 1:8
    [mdPrm_ctx(tt,:),~,sdPrm_ctx(tt,:)] = meanstdsem(cell2mat(dPrmTrjC(~[dPrmRez(:).isStr]',tt))); 
end
% plot pure load
figure; hold on; 
%plot(mean(mdPrm_ctx([4,8],:))./max(mean(mdPrm_ctx([4,8],:))))
plot(mdPrm_ctx(8,:)./max(mdPrm_ctx(8,:)))
plot(mdPrm_ctx(4,:)./max(mdPrm_ctx(4,:)))

% plot pure direction 
%plot(mean(mdPrm([2,6],:))./max(mean(mdPrm([2,6],:))))
plot(mean(mdPrm_ctx([2,6],:)))

% plot pure mixed 
%plot(mean(mdPrm([1,3,5,7],:))./max(mean(mdPrm([1,3,5,7],:))))
plot(mean(mdPrm_ctx([1,3,5,7],:)))
axis tight
set(gca,'tickDir','out')
print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData/collectFigure','dPrm_neuronType_Ctx'),'-dpdf','-painters')

