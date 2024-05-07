function imgFig= imagescWithTimeInfo(imgMatTrialByTime, varargin)

% Define the size of the matrix
[num_trials, num_time_bins] = size(imgMatTrialByTime);
y = 1:num_trials;
    
if nargin == 2 && length(varargin{1})==num_time_bins
    x = varargin{1}; 
else
    x = 1:num_time_bins; % Assuming your time bins are evenly spaced
end



% Create vectors for x and y coordinates


% Create the grid of x and y coordinates using meshgrid
%[X, Y] = meshgrid(x, y);

% Plot the matrix using imagesc
imgFig = figure; 
imagesc(x, y, imgMatTrialByTime);
colorbar; % Add colorbar for visualization
xlabel('Time Bins');
ylabel('Trials');

% Optionally, set x and y tick labels
% You can use imgMatTrialByTime_bins for x tick labels if needed
% set(gca, 'XTick', 1:numel(imgMatTrialByTime_bins), 'XTickLabel', imgMatTrialByTime_bins);
% set(gca, 'YTick', 1:num_trials, 'YTickLabel', 1:num_trials);

end