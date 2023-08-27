function [anova_dat, anova_label] = organize_outcell(outcell)

anova_dat = []; 
anova_label = {}; 
for jj = 1:size(outcell, 2)
    % get data
    dat_col = cell2mat(outcell(:, jj)); 
    dat_col = dat_col(~isnan(dat_col));
    anova_dat = [anova_dat; dat_col];
    
    % get label
    label_col = cell(length(dat_col), 1); 
    label_col(:) = {num2str(jj)}; 
    anova_label = [anova_label; label_col]; 
end

end