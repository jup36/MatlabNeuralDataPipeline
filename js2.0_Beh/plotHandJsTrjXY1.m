%% Helper functions
function [handXYuptoP1,jXY1] = plotHandJsTrjXY1(handXY,handXY1med,jsXY,jsXY1,rwdTrialI,colorTheme,figSavePath,figSaveName)

handXYuptoP1 = cell(size(handXY,3),1); 

ntb = size(handXY,2); % the # of timeBins
ntr = size(handXY,3); % the # of trials 

[cTheme] = TNC_CreateRBColormapJP(ntr*2,colorTheme); % color to assign across trials
%c = cTheme(max(1,ntr-floor(ntr/2)):max(1,ntr-floor(ntr/2))+ntr-1,:); % pick colors from the middle ones
c = cTheme(max(end-ntr+1-5,1):end-5,:); % pick colors from the middle ones

handXY = -handXY; 
handXY1med = -handXY1med;
jsXY = -jsXY;
jsXY1 = -jsXY1; 

%medHandXY1 = nanmedian(handXY(:,1,:),3); % median starting point (reference point)

hXY = handXY-repmat(handXY1med,[1,size(handXY,2),size(handXY,3)]); % normalized by subtracting the median initial hand position 
jXY = jsXY-repmat(handXY1med,[1,size(jsXY,2),size(jsXY,3)]); % normalized by subtracting the median initial hand position ; 
jXY1 = jsXY1-handXY1med; 

% get trajectories by the proximity to the joystick target position
[~,minDjXY1h] = min(sum((hXY-repmat(jXY1,[1,size(hXY,2),size(hXY,3)])).^2),[],2);
p1h = squeeze(minDjXY1h); % nth time bin that is most close to the initial joystick position (putative pullStart point)

[~,minDjXY1j] = min(sum((jXY-repmat(jXY1,[1,size(jXY,2),size(jXY,3)])).^2),[],2);
jp1 = squeeze(minDjXY1j); % nth time bin that is most close to the initial joystick position (putative pullStart point)

figure; hold on;
%scatter(jXY(1), jXY(2), 100, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0 0 0]) % draw joystick 
fstValTr = find(~isempty(hXY));
for j = 1:size(hXY,3) % trials to draw (just draw them all)
    if rwdTrialI(j)
        x = hXY(1,1:p1h(j),j);
        y = hXY(2,1:p1h(j),j);
    else
        x = hXY(1,:,j);
        y = hXY(2,:,j);
    end
    handXYuptoP1{j,1} = [x;y];
    % draw joystick initial position
    if j == fstValTr
        scatter(jXY1(1), jXY1(2), 200, 'MarkerEdgeColor', 'none','MarkerFaceColor','k','MarkerFaceAlpha',.7) % draw starting point hTrj
    end
    % draw initial and endpoint hand positions
    scatter(x(1), y(1), 50, 'MarkerEdgeColor', 'none','MarkerFaceColor',c(j,:),'MarkerFaceAlpha',.4) % draw starting point hTrj
    scatter(x(end), y(end), 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor',c(j,:),'MarkerFaceAlpha',.7) % draw last point hTrj
    plot(x,y,'color',c(j,:),'lineWidth',2)
    %patch([x nan],[y nan],[1:length(y) nan],'FaceColor','none','EdgeColor','interp','lineWidth',2)
    %colormap(c)
    
    % FOR JOYSTICK ALSO JUST DRAW FROM THE CLOSEST POINT FROM THE KNOWN jXY
    if j<=size(jXY,3)
        jx = jXY(1,jp1(j):end,j);
        jy = jXY(2,jp1(j):end,j);
    
        plot([jx(1) jx(end)], [jy(1) jy(end)],'k','lineWidth',2)
    end
end
xlim([-12 6])
ylim([-5 15])
pbaspect([1 1 1])
set(gca,'tickDir','out')
colormap(c); colorbar

print(fullfile(figSavePath,figSaveName),'-dpdf','-bestfit','-painters')
hold off;
end