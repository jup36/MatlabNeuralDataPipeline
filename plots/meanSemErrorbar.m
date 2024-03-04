function fig = meanSemErrorbar(mean_values, sem_values)

fig = figure; 
% Assuming mean_values and sem_values are already defined as 1x20 vectors

% Generate a vector for the x-axis
x_values = 1:length(mean_values); 

% Plot filled circles for mean values
scatter(x_values, mean_values, 'filled', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');

% Hold on to the current figure
hold on;

% Add error bars for standard error of the mean
errorbar(x_values, mean_values, sem_values, 'LineStyle', 'none', 'Color', 'r', 'CapSize', 10);

% Label the axes
xlabel('Dims');
ylabel('Mean Value');

% Add a title and a legend
%title('Mean Values with SEM');
%legend('Mean values', 'SEM', 'Location', 'best');

% Hold off to finish the plotting
hold off;

end