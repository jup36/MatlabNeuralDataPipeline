 filePaths = {
    '/Volumes/Extreme SSD/js2p0/WR38_052119/Matfiles', ... % Dual recording without Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR38_052419/Matfiles', ... % Corticostriatal recording M1 silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR39_100219/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_081919/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    '/Volumes/Extreme SSD/js2p0/WR40_082019/Matfiles', ... % Dual recording with contra Cg silencing (checked)
    };   

for f = 1:length(filePaths) 
    load(fullfile(filePaths{f}, 'magnitude_inhibition_breakThrough_gpfa'), 'rez')
    magBtC{f, 1} = rez.magBT; 
    magInhC{f, 1} = rez.magInh; 
    magReachC{f, 1} = rez.magReach; 
    sigInhIC{f, 1} = rez.sigInhI; 
    fprintf('Completed preprocessing file #%d\n', f)
end

magInhAll = cell2mat(magInhC);
magBtAll = cell2mat(magBtC);
magReachAll = cell2mat(magReachC); 
sigInhAll = cell2mat(sigInhIC); 

scatter_row_by_row_fillLogic([magInhAll, magBtAll, magReachAll], [-1.5 2.5], sigInhAll); 
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'cortical_inactivation_breakThrough_normalReach_allUnits'), '-dpdf', '-vector', '-bestfit')

scatter_row_by_row_fillLogic([magInhAll(sigInhAll), magBtAll(sigInhAll), magReachAll(sigInhAll)], [-1.5 2.5], ones(sum(sigInhAll), 1)); 
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'cortical_inactivation_breakThrough_normalReach_sigInhUnits'), '-dpdf', '-vector', '-bestfit')

%scatter_row_by_row_fillLogic([magInhC{f, 1}, magBtC{f, 1}, magReachC{f, 1}], [-1 3], sigInhIC{f, 1}); 

% Generate scatter plot with pastel colors and unity line
% Create a colormap with 100 pastel colors
cmap = pastelColormap(100);

% Assign colors to each point
numPoints = length(magBtAll); % Assuming magBtAll and magReachAll have the same length
colors = cmap(mod(0:numPoints-1, 100) + 1, :);

% Plot unity line (dotted)
figure; hold on;
maxVal = max([magBtAll; magReachAll]); % Find max value for plotting the unity line
minVal = min([magBtAll; magReachAll]); % Find min value for plotting the unity line
plot([minVal maxVal], [minVal maxVal], 'k--', 'LineWidth', 1.5); % Unity line

% Scatter plot with pastel colors
scatter(magBtAll, magReachAll, 40, colors, 'filled', 'MarkerFaceAlpha', 0.5);

% Make x and y axes equal
pbaspect([1 1 1])

% Add labels and grid
xlabel('Magnitude BT All');
ylabel('Magnitude Reach All');
set(gca, 'TickDir', 'out')
grid on;
title('Scatter Plot with Unity Line and Pastel Colors');
print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure', 'cortical_breakThrough_corticalReach_allUnits'), '-dpdf', '-vector', '-bestfit')


%% run repeated measures of ANOVA
%ibtr = [magInhAll, magBtAll, magReachAll]; 
ibtr = [magInhAll(sigInhAll), magBtAll(sigInhAll), magReachAll(sigInhAll)]; 


% artificial grouping to get by ranova
group_var = cell(size(ibtr, 1), 1); 
[group_var{1:ceil(size(ibtr, 1)/2)}] = deal('m1'); 
[group_var{ceil(size(ibtr, 1)/2)+1:end}] = deal('m2'); 

% variable names 
var_names = cell(1, 1+size(ibtr, 2)); 
for jj = 1:length(var_names)
    if jj == 1
        var_names{jj} = 'G'; 
    else
        var_names{jj} = strcat('B', num2str(jj-1)); 
    end
end

% build table
rm_tab = cell2table([group_var, num2cell(ibtr)], 'VariableNames', var_names); 

Block = table((1:size(ibtr,2))','VariableNames',{'Block'}); % Block
% Fit repetitive model to data
rm_design = ['B1-', strcat('B', num2str(size(ibtr,2))), ' ', '~', ' ', 'G'];  % rm design string

rm = fitrm(rm_tab, rm_design, 'WithinDesign', Block); % repeated measures model fit for RT data
ranova_rez = ranova(rm); % repeated measures ANOVA on RT (within- and within*group interaction)
multiCompBlock = multcompare(rm,'Block');

rez.rm = rm; 
rez.ranova_rez = ranova_rez; 
rez.multiCompBlock = multiCompBlock;  
rez.rm_tab = rm_tab; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create pastel colormap
function cmap = pastelColormap(nColors)
    baseColors = lines(nColors); % Start with MATLAB's 'lines' colormap
    cmap = 0.8 * baseColors + 0.2; % Blend with white to make pastel
end
