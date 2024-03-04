function rez = blockByblock_rm_simple(sbb0)
% sbb0: session-by-block data 
    
%% run repeated measures of ANOVA
% artificial grouping to get by ranova
group_var = cell(size(sbb0, 1), 1); 
[group_var{1:ceil(size(sbb0, 1)/2)}] = deal('m1'); 
[group_var{ceil(size(sbb0, 1)/2)+1:end}] = deal('m2'); 

% variable names 
var_names = cell(1, 1+size(sbb0, 2)); 
for jj = 1:length(var_names)
    if jj == 1
        var_names{jj} = 'G'; 
    else
        var_names{jj} = strcat('B', num2str(jj-1)); 
    end
end

% build table
rm_tab = cell2table([group_var, num2cell(sbb0)], 'VariableNames', var_names); 

Block = table((1:size(sbb0,2))','VariableNames',{'Block'}); % Block
% Fit repetitive model to data
rm_design = ['B1-', strcat('B', num2str(size(sbb0,2))), ' ', '~', ' ', 'G'];  % rm design string

rm = fitrm(rm_tab, rm_design, 'WithinDesign', Block); % repeated measures model fit for RT data
ranova_rez = ranova(rm); % repeated measures ANOVA on RT (within- and within*group interaction)
multiCompBlock = multcompare(rm,'Block');

rez.rm = rm; 
rez.ranova_rez = ranova_rez; 
rez.multiCompBlock = multiCompBlock;  
rez.rm_tab = rm_tab; 

end