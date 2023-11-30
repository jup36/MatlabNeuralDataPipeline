function histo = lickDist(tbytDat)
% This function examines the distribution of lick intervals
upperlim = 1; % the upper limit of lick intervals to be included 

lickIntervals = diff(cell2mat(cellfun(@(a) a', {tbytDat(:).Lick}, 'un', 0))); 
lickIntervals = lickIntervals(lickIntervals<upperlim); 

hold on; 
histo = histogram(lickIntervals, 50); 

cutoff = 0.4; 

title('Lick intervals')
xline(cutoff, 'r:', 'LineWidth', 1.5)
set(gca, 'TickDir', 'out')
xlabel('Interval Time (s)')
ylabel('Counts')

hold off; 

% print(fullfile('/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823/Figure', 'lick_interval_distribution_DA008_101823'), '-dpdf', '-vector')

end