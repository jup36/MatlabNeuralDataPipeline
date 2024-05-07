function ranges = calculateSymmetricPercentileRanges(data, step)
    % Calculate symmetric percentile ranges centered at zero, extending outward.
    %
    % Parameters:
    %   data - An array containing the data points.
    %   step - The step size for the percentiles in terms of percentage (e.g., 10 for 10%, 20% ...).
    %
    % Returns:
    %   ranges - A matrix where each row represents the [lower_bound, upper_bound] of the percentile range,
    %            including ranges from negative to positive.

    if nargin < 2
        step = 10; % Default step size is 10 if not specified
    end

    % Calculate the absolute values of percentiles from 50 to 100 in step increments
    upperPercentiles = 50:step:100;
    percentileValues = prctile(abs(data), upperPercentiles);

    % Calculate symmetric ranges
    numRanges = numel(percentileValues);
    ranges = zeros(2 * numRanges, 2);

    % Fill negative to positive ranges
    for i = 1:numRanges
        lowerBound = -percentileValues(numRanges - i + 1);
        upperBound = percentileValues(numRanges - i + 1);

        if i == 1
            ranges(2 * i - 1, :) = [lowerBound, 0];
            ranges(2 * i, :) = [0, upperBound];
        else
            ranges(2 * i - 1, :) = [lowerBound, -percentileValues(numRanges - i)];
            ranges(2 * i, :) = [percentileValues(numRanges - i), upperBound];
        end
    end

    % Sort the ranges for presentation
    ranges = sortrows(ranges, 1);
end
