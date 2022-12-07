%% Helper functions
function [val_trials] = plot_handJsTrjXY(filePath, fileName, varargin)

%filePath = fullfile('D:\Junchol_Data\JS2p0\WR40_081919\Matfiles','js2p0_tbytSpkHandJsTrjBin_WR40_081919.mat'); 
load(fullfile(filePath, fileName), 'ss', 'jkvt')

p = parse_input_plot_handJsTrjXY(filePath, fileName, varargin); 
%p = parse_input_plot_handJsTrjXY(filePath, {'draw_blocks', [1, 3], 'draw_trials', [30], ... 
%    'rewarded_only', true, 'draw_only_to_pullstart', true, 'draw_trials_per_block', 10}); 

%% define trials to draw
% pick by blocks
trI = cell2mat(cellfun(@(a) any(p.Results.draw_blocks(:)==a), {ss(:).blNumber}, 'un', 0));
% pick by trials
trI(p.Results.draw_trials) = true; 

% pick by rewarded or not
if p.Results.rewarded_only==true
   rewardI = [jkvt(:).rewarded]; 
   trI = trI & rewardI; 
end

% final trials to draw
val_trials = find(trI);
val_blocks = [jkvt(val_trials).blNumber]; 

% unique blocks 
blocks = unique([jkvt(val_trials).blNumber]); 

% get joystick trajectories
jXY = {ss(val_trials).jXY_n_r}; 

% get initial joystick positions
jXY1 = {ss(val_trials).jXY1_n_r}; 

% get hand trajectories 
if p.Results.draw_only_to_pullstart == true
    hXY = {ss(val_trials).hXY_to_p1}; % hand trajectories upto the pull start point
else
    hXY = {ss(val_trials).hXY_n_r};
end

%[cTheme] = TNC_CreateRBColormapJP(ntr*2, 'wblue'); % color to assign across trials

%% draw
figure; hold on;
% draw joystick initial positions
for b = blocks % 1:length(blocks)
    jXY1_block = median(cell2mat(jXY1(b==val_blocks)), 2); 
    scatter(jXY1_block(1), jXY1_block(2), 200, 'MarkerEdgeColor', 'none','MarkerFaceColor','k','MarkerFaceAlpha',.7) % draw starting point hTrj

    hXY_block = hXY(val_blocks==b); 
    
    % draw initial and endpoint hand positions
    if ismember(b, [1,2,5,6]) % left blocks
       c = [0, 0, 1]; % target on the LEFT
       if ismember(b, [1,5]) % left/light blocks
          width = 1;
       elseif ismember(b, [2,6]) % left/heavy blocks
          width = 2;
       end
    else % right blocks
       c = [1, 0, 0]; % target on the RIGHT
       if ismember(b, [3,7]) % right/light blocks
           width = 1; 
       elseif ismember(b, [4,8]) % right/heavy blocks
           width = 2; 
       end
    end
    
    % select hand trajectories to draw
    if ~isempty(p.Results.draw_trials_per_block) 
       if sum(val_blocks==b)> p.Results.draw_trials_per_block  
          block_dist_to_jXY1 = {ss(val_trials(val_blocks==b)).min_hXY_dist_to_jXY1}; %[ss(val_trials(val_blocks==b)).min_hXY_dist_to_jXY1];
          block_dist_to_jXY1_NaNs = cell2mat(cellfun(@isempty, block_dist_to_jXY1, 'un', 0)); 
          [block_dist_to_jXY1{block_dist_to_jXY1_NaNs}] = deal(NaN);
          [~,distI] = sort(cell2mat(block_dist_to_jXY1)); 
          hXY_b = hXY_block(distI(1:min(p.Results.draw_trials_per_block, sum(~block_dist_to_jXY1_NaNs)))); 
       elseif sum(val_blocks==b) <= p.Results.draw_trials_per_block  
          hXY_b = hXY_block; 
       end
    else % if there is no limit on trials to draw
        hXY_b = hXY_block; 
    end
    
    for j = 1:length(hXY_b) % draw trial by trial 
        if ~isempty(hXY_b{j})
            trj = hXY_b{j}; 
            scatter(trj(1,1), trj(2,1), 50, 'MarkerEdgeColor', 'none','MarkerFaceColor',c,'MarkerFaceAlpha',.4) % draw starting point hTrj
            scatter(trj(1,end), trj(2,end), 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor',c,'MarkerFaceAlpha',.7) % draw last point hTrj
            plot(trj(1,:), trj(2,:), 'color', c, 'lineWidth', width)
        end
    end
end

hold off;
mName = fileName(strfind(fileName,'WR'):end); 
print(fullfile(filePath, 'Figure', strcat('plot_handJsTrjXY_', mName)), '-dpdf', '-bestfit', '-painters')    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%
function p = parse_input_plot_handJsTrjXY(filePath, fileName, vargs)
% parse input, and extract name-value pairs
default_draw_blocks = []; %
default_draw_trials = []; %
default_numb_to_draw = 10;
default_rewarded_only = true;
default_draw_only_to_pullstart = false;
default_draw_trials_per_block = [];

p = inputParser; % create parser object
addRequired(p,'filePath');
addRequired(p,'fileName'); 
addParameter(p,'draw_blocks', default_draw_blocks);
addParameter(p,'draw_trials', default_draw_trials);
addParameter(p,'numb_to_draw', default_numb_to_draw);
addParameter(p,'rewarded_only', default_rewarded_only);
addParameter(p,'draw_only_to_pullstart', default_draw_only_to_pullstart)
addParameter(p,'draw_trials_per_block', default_draw_trials_per_block)

parse(p, filePath, fileName, vargs{:})
end

