function runRepeatedMeasuresANOVA(varargin)
    % Check if there are at least two conditions
    if nargin < 2
        error('At least two conditions are required for repeated measures ANOVA.');
    end
    
    % Number of conditions
    numConditions = nargin;
    
    % Combine all conditions into a table
    conditionData = table(varargin{:}, 'VariableNames', strcat('Condition', string(1:numConditions)));
    
    % Create a within-subject design table
    conditions = (1:numConditions)';
    withinDesign = table(conditions, 'VariableNames', {'Conditions'});
    
    % Fit the repeated measures model
    rm = fitrm(conditionData, strcat('Condition1-Condition', num2str(numConditions), ' ~ 1'), 'WithinDesign', withinDesign);
    
    % Run the repeated measures ANOVA
    ranovaResults = ranova(rm);
    
    % Display the ANOVA results
    disp('Repeated Measures ANOVA Results:');
    disp(ranovaResults);
    
    % Perform multiple comparisons
    multComp = multcompare(rm, 'Conditions', 'ComparisonType', 'tukey-kramer');
    
    % Display the multiple comparison results
    disp('Multiple Comparisons Results (Tukey HSD):');
    disp(multComp);
end
