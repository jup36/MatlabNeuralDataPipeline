%% Helper functions
function [handXYuptoP1,jXY1] = plotHandJsTrjXY1_clockwiseRotation(filePath, figSavePath, figSaveName, theta_in_deg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filePath = 'D:\Junchol_Data\JS2p0\WR40_081919\Matfiles','js2p0_tbytSpkHandJsTrjBin_WR40_081919.mat'; 
load(fullfile(filePath), 'ss', 'jkvt')
valTrs = cell2mat(cellfun(@(a) ~isempty(a), {ss(:).hTrjB}, 'un', 0));
rwdTrI = [jkvt(:).rewarded]; 


valrwdTrs = find(valTrs & rwdTrI); % valid rewarded trials 

% get hand and joystick XY trajectories of select trials
hXY_C = cellfun(@(a) -a(1:2,:), {ss(valTrs & rwdTrI).hTrjB}, 'un', 0); 
hXY1med = median(cell2mat({ss(valTrs & rwdTrI).hInitPos}),2); 
hXY1med = -hXY1med(1:2,:); 

jXY_C = cellfun(@(a) -a(1:2,:), {ss(valTrs & rwdTrI).jTrjB}, 'un', 0);
jXY1med = median(cell2mat({ss(valTrs & rwdTrI).jsXYbot}),2);
jXY1med = -jXY1med(1:2,:); 

% get rotation matrix
%theta_in_deg = 15; % in degree
theta = deg2rad(theta_in_deg);   
cwRot2d = [cos(theta), sin(theta); -sin(theta), cos(theta)]; 
% cwRot3d = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1]; % use this for 3D trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessing - baseline (hand initial position) subtraction and 2D rotation                                                  
baseSub_rotateXY_hTrj_f = @(x) ((x-hXY1med)'*cwRot2d)'; % function handle for hand trajectory baseline subtraction and rotation

hXY_n_r = cellfun(@(a) baseSub_rotateXY_hTrj_f(a), hXY_C, 'un', 0); % baseline subtraction and rotation
hXY1_n_r = ((hXY1med-hXY1med)'*cwRot2d)'; 

jXY_n_r = cellfun(@(a) baseSub_rotateXY_hTrj_f(a), jXY_C, 'un', 0); % baseline subtraction and rotation
jXY1_n_r = ((jXY1med-hXY1med)'*cwRot2d)'; 

% preprocessing - identify the point when hand was closest to the joystick initial position
dist_to_jXY1_f = @(x) sqrt(sum((x-repmat(jXY1_n_r, [1, size(x,2)])).^2)); 
minI_f = @(x) find(x==min(x)); 

hXY_dist_to_j1 = cellfun(@(h) dist_to_jXY1_f(h), hXY_n_r, 'un', 0); 

p1_bin = cellfun(@(h) minI_f(h), hXY_dist_to_j1, 'un', 0); 
hXY_to_p1 = cellfun(@(a,b) a(:,1:b), hXY_n_r, p1_bin , 'un', 0); 



cellfun()


min_dist_to_jXY1_f = @(x) min(sqrt(sum((x-repmat(jXY1_n_r, [1, size(x,2)])).^2)), [], 2);  
hXY_min_dist_to_j1 = cell2mat(cellfun(@(h) min_dist_to_jXY1_f(h), hXY_n_r, 'un', 0)); 











ntb = size(handXY,2); % the # of timeBins
ntr = size(handXY,3); % the # of trials 

[cTheme] = TNC_CreateRBColormapJP(ntr*2,colorTheme); % color to assign across trials
%c = cTheme(max(1,ntr-floor(ntr/2)):max(1,ntr-floor(ntr/2))+ntr-1,:); % pick colors from the middle ones
c = cTheme(max(end-ntr+1-5,1):end-5,:); % pick colors from the middle ones






% get trajectories by the proximity to the joystick target position
[~,minDjXY1h] = min(sum((hXY-repmat(jXY1,[1,size(hXY,2),size(hXY,3)])).^2),[],2);
p1h = squeeze(minDjXY1h); % nth time bin that is most close to the initial joystick position (putative pullStart point)

[~,minDjXY1j] = min(sum((jXY-repmat(jXY1,[1,size(jXY,2),size(jXY,3)])).^2),[],2);
jp1 = squeeze(minDjXY1j); % nth time bin that is most close to the initial joystick position (putative pullStart point)

figure; hold on;
%scatter(jXY(1), jXY(2), 100, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor',[0 0 0]) % draw joystick 
fstValTr = find(~isempty(hXY));
for j = 1:size(hXY,3) % trials to draw (just draw them all)
    if rwdTrI(j)
        x = hXY(1,1:p1h(j),j);
        y = hXY(2,1:p1h(j),j);
    else
        x = hXY(1,:,j);
        y = hXY(2,:,j);
    end
    handXYuptoP1{j,1} = [x;y];
    handXYuptoP1{j,2} = cwRot2d*[x;y];
    
    xR = handXYuptoP1{j,2}(1,:); % rotated x
    yR = handXYuptoP1{j,2}(2,:); % rotated y
    
    % draw joystick initial position
    if j == fstValTr
        scatter(jXY1R(1), jXY1R(2), 200, 'MarkerEdgeColor', 'none','MarkerFaceColor','k','MarkerFaceAlpha',.7) % draw starting point hTrj
    end
    % draw initial and endpoint hand positions
    scatter(xR(1), yR(1), 50, 'MarkerEdgeColor', 'none','MarkerFaceColor',c(j,:),'MarkerFaceAlpha',.4) % draw starting point hTrj
    scatter(xR(end), yR(end), 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor',c(j,:),'MarkerFaceAlpha',.7) % draw last point hTrj
    plot(xR,yR,'color',c(j,:),'lineWidth',2)
    %patch([x nan],[y nan],[1:length(y) nan],'FaceColor','none','EdgeColor','interp','lineWidth',2)
    %colormap(c)
    
    % FOR JOYSTICK ALSO JUST DRAW FROM THE CLOSEST POINT FROM THE KNOWN jXY
    if j<=size(jXY,3)
        jx = jXY(1,jp1(j):end,j);
        jy = jXY(2,jp1(j):end,j);
        jxy = [jx;jy];
        jxyr = cwRot2d*jxy; 
        
        plot([jxyr(1,1) jxyr(1,end)], [jxyr(2,1) jxyr(2,end)],'k','lineWidth',2)
    end
end
xlim([-12 6])
xticks(-15:5:10)
ylim([-5 15])
yticks(-5:5:15)
pbaspect([1 1 1])
set(gca,'tickDir','out')
colormap(c); colorbar

print(fullfile(figSavePath,figSaveName),'-dpdf','-bestfit','-painters')
hold off;
end