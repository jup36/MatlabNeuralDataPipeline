% Force trajectory aligned to the max pull (negative) force point  
function [forceTrjCollect] = plotForceAlignToMaxTorque(forceTrjC, figSavePath, figSaveName, nBackFromMaxForce, colorScheme)
% forceTrjC = {ss(trials1).forceB};
% nBackFromMaxForce = 11; 
% colorScheme = 'wblue'; 
% figSaveName = 'forceTrajectories_b2_first10'; 
nTr = length(forceTrjC); 
forceTrjCollect = nan(nTr,nBackFromMaxForce); 
figure; hold on; 
[cTheme] = TNC_CreateRBColormapJP(nTr*2,colorScheme); % color to assign across trials
c = cTheme(end-nTr+1:end,:); % pick colors from the middle ones
for j = 1:nTr
    if ~isempty(forceTrjC{j})
        [minV,minI]=min(forceTrjC{j}); 
        plot(forceTrjC{j}(max(1,minI-nBackFromMaxForce+1):minI),'color',c(j,:))
        scatter(nBackFromMaxForce, minV, 100, 'MarkerEdgeColor', 'none','MarkerFaceColor',c(j,:),'MarkerFaceAlpha',.8) % draw starting point hTrj
        forceTrjCollect(j,:)=forceTrjC{j}(max(1,minI-nBackFromMaxForce+1):minI); 
        %forceTrjCollect(j,end-length(forceTrjC{j}(max(1,minI-nBackFromMaxForce+1):minI))+1:end)=forceTrjC{j}(max(1,minI-nBackFromMaxForce+1):minI); 
    end
end
plot(nanmean(forceTrjCollect),'color',cTheme(end,:),'lineWidth',3)
xlim([1 size(forceTrjCollect,2)+1])
ylim([-120 0])
hold off;
print(fullfile(figSavePath,figSaveName),'-dpdf','-bestfit','-painters')
end 