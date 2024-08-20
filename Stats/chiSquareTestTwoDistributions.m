function [pValue, chi2stat, df] = chiSquareTestTwoDistributions(distribution1, distribution2, groupIndex)

assert(length(distribution1) == length(distribution2), 'Distributions must be of equal length');

% Create the contingency table
contingencyTable = [sum(distribution1(groupIndex)), sum(distribution1) - sum(distribution1(groupIndex)); 
                    sum(distribution2(groupIndex)), sum(distribution2) - sum(distribution2(groupIndex))];

% Perform Chi-squared test
[chi2stat, pValue, df] = chi2cont(contingencyTable);

% Display the result
%fprintf('Chi-squared test for group %d: p-value = %.4f, chi2 = %.4f, df = %d\n', groupIndex, pValue, chi2stat, df);

end

function [chi2stat, pValue, df] = chi2cont(contingencyTable)
    % Get the row and column sums
    rowSums = sum(contingencyTable, 2);
    colSums = sum(contingencyTable, 1);
    total = sum(contingencyTable(:));

    % Compute the expected frequencies
    expected = rowSums * colSums / total;

    % Compute the chi-squared statistic
    chi2stat = sum((contingencyTable - expected).^2 ./ expected, 'all');

    % Compute the degrees of freedom
    df = (size(contingencyTable, 1) - 1) * (size(contingencyTable, 2) - 1);

    % Compute the p-value
    pValue = 1 - chi2cdf(chi2stat, df);
end
