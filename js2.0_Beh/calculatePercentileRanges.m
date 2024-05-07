function ranges = calculatePercentileRanges(data, step)
    % Calculate percentile ranges for given step sizes in a dataset.
    %
    % Parameters:
    %   data - An array containing the data points.
    %   step - The step size for the percentiles (e.g., 10 for 10%, 20% ...).
    %
    % Returns:
    %   ranges - A matrix where each row represents the [lower_bound, upper_bound] of the percentile range.

    if nargin < 2
        step = 10; % Default step size is 10 if not specified
    end
    
    percentiles = 0:step:100; % Create a vector of percentiles from 0 to 100
    values = prctile(data, percentiles); % Calculate the percentiles
    
    % Initialize the ranges matrix
    numRanges = length(values) - 1;
    ranges = zeros(numRanges, 2);
    
    for i = 1:numRanges
        ranges(i, 1) = values(i);
        ranges(i, 2) = values(i + 1);
    end
end