function plot_within_reach_vs_pull(reach_cell, pull_cell, region, ylimits)
regionI = cell2mat(cellfun(@(a) strcmp(a, region), reach_cell(1,:), 'un', 0));

reach_region = reach_cell(2:end, regionI); 
pull_region = pull_cell(2:end, regionI); 

collect_C = [reach_region, pull_region]; 
empties = cell2mat(cellfun(@isempty, collect_C, 'un', 0));

[collect_C{empties}] = deal(NaN); 
collect_region = cell2mat(cellfun(@(a) a(end), collect_C, 'un', 0)); 

scatter_row_by_row(collect_region, ylimits)

end