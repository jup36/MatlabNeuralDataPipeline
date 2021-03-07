function plotHandAtInitReachPosition(initHcell, reachHcell, jXY1cell, figSavePath, saveNameCell)
% initHcell = stat.initHcell;
% reachHcell = stat.reachHcell;
% jXY1cell = {rez.jXY1};
lb = [1,2,5,6]; % left target position blocks
rb = [3,4,7,8]; % right target position blocks

for f = 1:size(initHcell,1) % file
    figure; hold on;
    for b = 1:size(initHcell,2) % block
        if ~isempty(initHcell{f,b,1})&&~isempty(reachHcell{f,b,1})
            if ismember(b,lb)
                c = 'b';
            else
                c = 'r';
            end
            tmpHInit = cell2mat(squeeze(initHcell(f,b,:))');
            tmpHReach = cell2mat(squeeze(reachHcell(f,b,:))');
            tmpJ = jXY1cell{f}{b};
            
            % draw joystick initial position
            scatter(tmpJ(1), tmpJ(2), 300, 'MarkerEdgeColor', 'none','MarkerFaceColor','k','MarkerFaceAlpha',.7) % draw starting point hTrj
            % draw initial and endpoint hand positions
            scatter(tmpHInit(1,:), tmpHInit(2,:), 50, 'MarkerEdgeColor', c,'MarkerFaceColor','none','MarkerEdgeAlpha',.4) % draw starting point hTrj
            scatter(tmpHReach(1,:), tmpHReach(2,:), 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor',c,'MarkerFaceAlpha',.7) % draw last point hTrj
        end
    end
    hold off;
    pbaspect([1 1 1])
    set(gca,'tickDir','out')
    reachHall = cell2mat(reshape(squeeze(reachHcell(f,:,:)),1,[])); 
    reachmed = nanmedian(reachHall'); 
    reachstd = nanstd(reachHall'); 
    ylim([-5 reachmed(2)+3*reachstd(2)])
    print(fullfile(figSavePath,strcat('handScatterAtInitReachJs_',saveNameCell{f})),'-dpdf','-painters','-bestfit')
    %close all
end
end