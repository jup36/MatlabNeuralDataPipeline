function one = onewayANOVA_multiCompare(dat)
%This function runs one-way anova and posthoc pairwise test with correction
% for multiple comparisons.
% Input: dat is n-by-p matrix (n: data points, p: groups).
% Output: one (output structure). Note that the multiple comparison result
% is contained in one.c, see below for interpretation of c.
%For example, suppose one row contains the following entries.
% 2.0000  5.0000  1.9442  8.2206  14.4971 0.0432
% These numbers indicate that the mean of group 2 minus the mean of group 5 is estimated to be 8.2206,
% and a 95% confidence interval for the true difference of the means is [1.9442, 14.4971].
% The p-value for the corresponding hypothesis test that the difference of the means of groups 2 and 5
% is significantly different from zero is 0.0432.

% One-way ANOVA
[one.p, one.tbl, one.stats] = anova1(dat);
[one.c, ~, ~, ~] = multcompare(one.stats, "CriticalValueType", "tukey-kramer");
end