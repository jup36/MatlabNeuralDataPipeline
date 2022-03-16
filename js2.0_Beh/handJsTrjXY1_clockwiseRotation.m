%% Helper functions
function handJsTrjXY1_clockwiseRotation(filePath, fileName, theta_in_deg)

% filePath = 'D:\Junchol_Data\JS2p0\WR40_081919\Matfiles';
% fileName = 'js2p0_tbytSpkHandJsTrjBin_WR40_081919';
% theta_in_deg = -18; 
load(fullfile(filePath, fileName), 'ss')
valTrs = cell2mat(cellfun(@(a) ~isempty(a), {ss(:).hTrjB}, 'un', 0));
valTrs = find(valTrs); % valid rewarded trials
%rwdTrI = [jkvt(:).rewarded];

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
end

save(fullfile(filePath, fileName), 'ss', '-append')

end