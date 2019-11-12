function plotRepresentativeTrjsPerReachPos(filePath)

cSeed{2,1} = [0 0 205]./255; % B
cSeed{1,1} = [100 149 237]./255;

cSeed{2,2} = [255 0 0]./255; % R
cSeed{1,2} = [255 193 193]./255;

cSeed{2,3} = [0 255 0]./255; % G
cSeed{1,3} = [152 251 152]./255;

cSeed{2,4} = [255 0 255]./255; % M
cSeed{1,4} = [255 225 255]./255;

numbTrialsToPlot = 10; 

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

jsTargetP = sort(unique([jkvt(:).reachP1]));
jsPullTrq = sort(unique([jkvt(:).pull_torque]));

%% sort/select trials to plot
rIdx = [jkvt(:).rewarded]; % reward index
hIdx = ~cellfun(@isempty, {jkvt(:).hTrjF}); % hand trajectory availability index

%figure; hold on;
% plot trajectories
for p = 1:length(jsTargetP)
    figure; hold on;
    
    for rp = 1:length(jsTargetP)
        rpTr = find(cellfun(@(a) a==jsTargetP(rp), {jkvt(:).reachP1}), 1, 'first');
        jsAtTarget = [jkvt(rpTr).jsTreachPosT jkvt(rpTr).jsTreachPosB]; % joystick at the reach target position
        jsX = jsAtTarget(1,:); % Js at the reach target position X
        jsY = jsAtTarget(2,:); % Js at the reach target position Y
        jsZ = jsAtTarget(3,:); % Js at the reach target position Z
        
        if p==rp
            plot3(jsX, jsY, jsZ, 'color', 'k', 'lineWidth', 3)
        else
            plot3(jsX, jsY, jsZ, 'color', [158 158 158]./255, 'lineWidth', 3)
        end
    end
    
    pIdx = [jkvt(:).reachP1]==jsTargetP(p); % position index
    cmod = mod(p,size(cSeed,2)); % pick color
    if cmod == 0
        cmod = size(cSeed,2);
    end
    [tmpJsC] = colorGradient(cSeed{1,cmod}, cSeed{2,cmod},length(jsPullTrq));
    
    for q = 1:length(jsPullTrq)
        qIdx = [jkvt(:).pull_torque]==jsPullTrq(q); % torque index
        c = tmpJsC(q,:);
        
        tmpTrs = find(pIdx & qIdx & rIdx & hIdx); % trials to plot individually
        [~,tmpTrsSrtI] = sort(cellfun(@(a) size(a,2), {jkvt(tmpTrs).hTrjF}));
        sortTrs = tmpTrs(tmpTrsSrtI);
        
        trsToPlot = sortTrs(1:min(numbTrialsToPlot,length(sortTrs)));
        
        for t = 1:length(trsToPlot)
            tt = trsToPlot(t);
            rStart = jkvt(tt).hTrjRstart(end);
            rStop = jkvt(tt).hTrjRstop(end);
            
            if ~isempty(rStart) && ~isempty(rStop) && rStart<rStop
                x = jkvt(tt).hTrjF(1,rStart:rStop);
                y = jkvt(tt).hTrjF(2,rStart:rStop);
                z = jkvt(tt).hTrjF(3,rStart:rStop);
                cr = 1:size(jkvt(tt).hTrjF(rStart:rStop),2); % generate a colormap;
                patch([x nan],[y nan],[z nan],[cr nan],'EdgeColor',c)
            end
        end
        clearvars t
    end
    xlabel('X(mm)')
    ylabel('Y(mm)')
    zlabel('Z(mm)')
    view(-142,41)
    hold off
    print(fullfile(filePath,'Figure',sprintf('representativeTrajOfReachPos%d',jsTargetP(p))),'-dpdf','-painters','-bestfit');
end
clearvars p 


end