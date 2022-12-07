function rez = blockByblock_rm(sbb0, file_path, varargin)
% sbb0: session-by-block behavioral data 

% preprocess the data matrix
if size(sbb0, 1) ~= length(file_path)
    error("Check the input data, they don't match!")
end

p = parse_input_blockByblock_rm(sbb0, file_path, varargin);

if p.Results.avg_m  % animal average
    for i = 1:length(file_path)
        wrI = strfind(file_path{i}, 'WR');
        m_name_c{i, 1} = file_path{i}(wrI:wrI+3);  
    end
    
    unique_m_name_c = unique(m_name_c); 
    sbb1 = nan(size(unique_m_name_c, 1), size(sbb0, 2)); 
    
    for j = 1:length(unique_m_name_c)
        m_I = cellfun(@(a) strcmpi(a, unique_m_name_c{j}), m_name_c); 
        sbb1(j, :) = nanmean(sbb0(m_I, :), 1);  
    end
else  % without animal average 
    sbb1 = sbb0; 
end

if p.Results.avg_b  % block average
   sbb2 = nan(size(sbb1, 1), length(p.Results.b_pairs)); 
   for j = 1:length(p.Results.b_pairs) 
       sbb2(:, j) = nanmean(sbb1(:, p.Results.b_pairs{j}), 2); 
   end
else  % without block average
    sbb2 = sbb1; 
end
    
%% run repeated measures of ANOVA
% artificial grouping to get by ranova
group_var = cell(size(sbb2, 1), 1); 
[group_var{1:ceil(size(sbb2, 1)/2)}] = deal('m1'); 
[group_var{ceil(size(sbb2, 1)/2)+1:end}] = deal('m2'); 

% variable names 
var_names = cell(1, 1+size(sbb2, 2)); 
for jj = 1:length(var_names)
    if jj == 1
        var_names{jj} = 'G'; 
    else
        var_names{jj} = strcat('B', num2str(jj-1)); 
    end
end

% build table
rm_tab = cell2table([group_var, num2cell(sbb2)], 'VariableNames', var_names); 

Block = table((1:size(sbb2,2))','VariableNames',{'Block'}); % Block
% Fit repetitive model to data
rm_design = ['B1-', strcat('B', num2str(size(sbb2,2))), ' ', '~', ' ', 'G'];  % rm design string

rm = fitrm(rm_tab, rm_design, 'WithinDesign', Block); % repeated measures model fit for RT data
ranova_rez = ranova(rm); % repeated measures ANOVA on RT (within- and within*group interaction)
multiCompBlock = multcompare(rm,'Block');

rez.rm = rm; 
rez.ranova_rez = ranova_rez; 
rez.multiCompBlock = multiCompBlock;  
rez.m_names = unique_m_name_c; 
rez.rm_tab = rm_tab; 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = parse_input_blockByblock_rm(sbb0, file_path, vargs)
        % sbb0: session-by-block behavioral data 
        default_avg_b = true;
        default_avg_m = true;
        default_b_pairs = {[1, 5], [2, 6], [3, 7], [4, 8]}; 
        
        p = inputParser; % create parser object
        addRequired(p, 'sbb0')
        addRequired(p, 'file_path')
        addParameter(p, 'avg_b', default_avg_b)
        addParameter(p, 'avg_m', default_avg_m)
        addParameter(p, 'b_pairs', default_b_pairs)
        
        parse(p, sbb0, file_path, vargs{:})
        
    end
end