%% Helper functions
function handJsTrjXY1_clockwiseRotation_tortuosity(filePath, theta_in_deg)

cd(filePath)
file = dir('**/*js2p0_tbytSpkHandJsTrjBin_*.mat');
load(fullfile(file.folder, file.name), 'ss')

valTrs = cell2mat(cellfun(@(a) ~isempty(a), {ss(:).hTrjB}, 'un', 0));
valTrs = find(valTrs); % valid rewarded trials

%valTrs_js = cell2mat(cellfun(@(a) ~isempty(a), {ss(:).jTrjB}, 'un', 0));
%valTrs_js = find(valTrs_js); % valid rewarded trials

% get hand and joystick XY trajectories of select trials
hXY_C = cellfun(@(a) -a(1:2,:), {ss(valTrs).hTrjB}, 'un', 0);
hXY1med = nanmedian(cell2mat({ss(valTrs).hInitPos}),2);
hXY1med = -hXY1med(1:2,:);

for j = unique([ss(:).blNumber]) 
    blockI = [ss(:).blNumber] == j; 
    js_init = [ss(blockI).jsXYbot];
    js_init_C = num2cell(repmat(js_init(:,1), 1, sum(blockI)),1); 
    [ss(blockI).jsXYbot] = js_init_C{:};  
end

for t = valTrs
    if isempty(ss(t).jTrjB)
        ss(t).jTrjB = ss(t).jsXYbot; 
    end   
end

jXY_C = cellfun(@(a) -a(1:2,:), {ss(valTrs).jTrjB}, 'un', 0);
jXY1_C = cellfun(@(a) -a(1:2,:), {ss(valTrs).jsXYbot}, 'un', 0);

% get rotation matrix
%theta_in_deg = 15; % in degree
theta = deg2rad(theta_in_deg);
cwRot2d = [cos(theta), sin(theta); -sin(theta), cos(theta)];
% cwRot3d = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1]; % use this for 3D trajectories

% preprocessing - baseline (hand initial position) subtraction and 2D rotation
baseSub_rotateXY_hTrj_f = @(x) ((x-hXY1med)'*cwRot2d)'; % function handle for hand trajectory baseline subtraction and rotation

hXY_n_r = cellfun(@(a) baseSub_rotateXY_hTrj_f(a), hXY_C, 'un', 0); % baseline subtraction and rotation

jXY_n_r = cellfun(@(a) baseSub_rotateXY_hTrj_f(a), jXY_C, 'un', 0); % baseline subtraction and rotation
jXY1_n_r = cellfun(@(a) baseSub_rotateXY_hTrj_f(a), jXY1_C, 'un', 0); % baseline subtraction and rotation

% preprocessing - identify the point when hand was closest to the joystick initial position
dist_to_jXY1_f = @(x, y) sqrt(sum((x-repmat(y, [1, size(x,2)])).^2));
minI_f = @(x) find(x==min(x));

hXY_dist_to_j1 = cellfun(@(h, j) dist_to_jXY1_f(h, j), hXY_n_r, jXY1_n_r, 'un', 0);

p1_bin = cellfun(@(h) minI_f(h), hXY_dist_to_j1, 'un', 0);
hXY_to_p1 = cellfun(@(a,b) a(:,1:b), hXY_n_r, p1_bin , 'un', 0);

for j = 1:length(valTrs)
    ss(valTrs(j)).hXY_n_r = hXY_n_r{j}; % recentered, rotated 2d hand trajectory
    ss(valTrs(j)).jXY_n_r = jXY_n_r{j}; % recentered, rotated 2d joystick trajectory
    ss(valTrs(j)).jXY1_n_r = jXY1_n_r{j}; % recentered, rotated 2d joystick initial position
    ss(valTrs(j)).hXY_to_p1 = hXY_to_p1{j}; % recentered, rotated 2d hand trajectory upto pull start
    ss(valTrs(j)).min_hXY_dist_to_jXY1 = min(hXY_dist_to_j1{j}); % minimum distance to the joystick initial position
    ss(valTrs(j)).rchAngDeg = computeReachAngle(hXY_n_r(j), jXY1_n_r{j}); 
    if ~isempty(ss(valTrs(j)).tPullStart) && ~isempty(hXY_to_p1{j})
        % get tortuosity = distance (path length) / displacement (endpoints)
        if size(hXY_to_p1{j},2)>50 %nBaseBins = 50
            reachTrj = hXY_to_p1{j}(:,max(1,50-3):end); % tightly (right before the reachStart) around the reach portion of the trajectory
        else
            reachTrj = hXY_to_p1{j}(:,1:end);
        end
        displace = sqrt(sum((reachTrj(:,end)-reachTrj(:,1)).^2)); % displacement: equivalently, sum(sqrt(sum(diff(shortTrj).^2,2)))
        distance = sum(sqrt(sum(diff(reachTrj).^2,2))); % actual distance
        ss(valTrs(j)).hXYtort = min(15,distance/displace); % tortuosity 
    end
end

save(fullfile(file.folder, file.name), 'ss', '-append')

end

function reachAngle = computeReachAngle( hTrjC, jsXYpos )
%Angle is calculated relative to the straight line projected from the
% initial hand position.  
    distToJs = cellfun(@(a) sum((a-repmat(jsXYpos,1,size(a,2))).^2,1), hTrjC, 'un', 0);  % compute the distance from the jsXY
    sign = [1,-1]; 
    for j = 1:length(distToJs)
        [~,tmpMinI] = min(distToJs{j});
        u=hTrjC{j}(:,tmpMinI)-hTrjC{j}(:,1); % reach vector
        %v=[hTrjC{j}(1,tmpMinI);hTrjC{j}(2,1)]-hTrjC{j}(:,1); % reference vector
        v=[hTrjC{j}(1,1);hTrjC{j}(2,tmpMinI)]-hTrjC{j}(:,1); % reference vector
        reachAngle(j) =  sign((u(1)>0)+1) * angleTwoVectors(u,v); % get the angle between reach and reference vectors      
    end
end
    