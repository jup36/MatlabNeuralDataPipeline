function plot3dHandTrajWithJsPos( filePath, trI)

% ensure that the jkvt structure is in the workspace
if exist('jkvt','var')
    
else
    jsKinFile = dir(fullfile(filePath,'jsTime1k_KinematicsTrajectories.mat'));
    
    if length(jsKinFile)==1
        jkvt = load(fullfile(filePath,'jsTime1k_KinematicsTrajectories.mat'),'jkvt'); % load jkvt
    else
        disp('Point to "jsTime1k_KinematicsTrajectories.mat" with the variable "jkvt"!')
        [jsKinFileSelect,jsKinPathSelect] = uigetfile(filePath);
        jkvt = load(fullfile(jsKinPathSelect,jsKinFileSelect),'jkvt'); % load jkvt
    end
    
    jkvt = jkvt.('jkvt'); 
end

% draw joystick target positions first 
jsTargetPos = unique([jkvt(:).reachP1]); 
jsC = cbrewer('seq', 'Greys', length(jsTargetPos)+3); 

figure; hold on; 
for rp = 1:length(jsTargetPos)
    rpTr = find(cellfun(@(a) a==jsTargetPos(rp), {jkvt(:).reachP1}), 1, 'first'); 
    jsAtTarget = [jkvt(rpTr).jsTreachPosT jkvt(rpTr).jsTreachPosB]; % joystick at the reach target position
    jsX = jsAtTarget(1,:); % Js at the reach target position X
    jsY = jsAtTarget(2,:); % Js at the reach target position Y
    jsZ = jsAtTarget(3,:); % Js at the reach target position Z
    
    if jkvt(rpTr).reachP1==jkvt(trI).reachP1 
        plot3(jsX, jsY, jsZ, 'color', 'k', 'lineWidth', 3)
    else
        plot3(jsX, jsY, jsZ, 'color', jsC(2,:), 'lineWidth', 3)
    end
end
    
% draw the hand trajectory
x = jkvt(trI).hTrjF(1,:); %T(jkvt(trI)).vUseFrameIdx);
y = jkvt(trI).hTrjF(2,:); %T(jkvt(trI)).vUseFrameIdx);
z = jkvt(trI).hTrjF(3,:); %T(jkvt(trI)).vUseFrameIdx);
c = 1:size(jkvt(trI).hTrjF,2); % generate a colormap;

patch([x nan],[y nan],[z nan],[c nan],'FaceColor','none','EdgeColor','interp')

xlabel('X(mm)')
ylabel('Y(mm)')
zlabel('Z(mm)')

colormap cool
colorbar
caxis([1 size(jkvt(trI).hTrjF,2)])
%print(fullfile(filePath,'Figure',sprintf('plot3dHandTrajWithJsPosTr#%d',trI)),'-dpdf','-painters','-bestfit');
hold off; 

end