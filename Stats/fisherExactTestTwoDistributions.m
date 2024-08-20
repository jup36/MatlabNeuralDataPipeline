function [pValue, stats] = fisherExactTestTwoDistributions(distribution1, distribution2, groupIndex)

assert(length(distribution1)==length(distribution2))

% Create the contingency table
contingencyTable = [sum(distribution1(groupIndex)), sum(distribution1) - sum(distribution1(groupIndex)); 
                    sum(distribution2(groupIndex)), sum(distribution2) - sum(distribution1(groupIndex))];

% Perform Fisher's exact test
[~, pValue, stats] = fishertest(contingencyTable);

% Display the result
%fprintf('Fisher''s exact test for group %d: p-value = %.4f\n', groupIndex, pValue);

end