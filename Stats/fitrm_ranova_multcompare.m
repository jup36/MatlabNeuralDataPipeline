% generate random data for the example
data = randn(168,3);
% define a variable that stores the treatment information
Treatment = cell(168, 1);
for i=1:25;Treatment{i} = 'control';end;
for i=26:52;Treatment{i} = 'treatment';end;
for i=53:78;Treatment{i} = 'control';end;
for i=79:98;Treatment{i} = 'treatment';end;
for i=99:137;Treatment{i} = 'control';end;
for i=138:168;Treatment{i} = 'treatment';end;
% define a variable that stores the environment type
Etype = cell(168, 1);
for i=1:52;Etype{i} = 'A';end;
for i=53:98;Etype{i} = 'B';end;
for i=99:168;Etype{i} = 'C';end;
% Store the data in a proper table format to do repeated measures analysis
all_table = table(Etype, Treatment, data(:, 1), data(:, 2), ...
    data(:, 3), ...
    'VariableNames', {'Etype', 'Treatment', 'pre', 'during', 'post'});
% Define the within subject parameter (pre, during, post)
Time = [0 1 2]; % 0: pre, 1: during, 2: post
% Fit repetitive model to data
rm = fitrm(all_table, 'pre-post ~ Treatment*Etype', ...
    'WithinModel', Time, 'WithinModel', 'separatemeans');
% Run repetitive measures ANOVA on model:
ranova(rm)
% make pairwise comparisons for the two-way interactions
multcompare(rm, 'Time', 'By', 'Etype')
multcompare(rm, 'Time', 'By', 'Treatment')
multcompare(rm, 'Treatment', 'By', 'Etype')
% but how can I make pairwise comparisons for the 3-way interaction?
% These do not work:
%multcompare(rm, 'Time', 'By', 'Treatment', 'Etype')
%multcompare(rm, 'Time', 'By', {'Treatment', 'Etype'})

% 1. Convert the between-subject factors to categorical.
etype_cat=categorical(Etype); treatment_cat=categorical(Treatment);

% 2. Create an interaction factor capturing each combination of levels % of Etype and Treatment (you can check with the function "catogories()")
interaction_cat=etype_cat.*treatment_cat;

% 3. Call fitrm with the modified between design.
t2 = table(interaction_cat, data(:, 1), data(:, 2), data(:, 3), 'VariableNames', {'interaction_etype_treat', 'pre', 'during', 'post'});   
rm2 = fitrm(t2, 'pre-post ~ interaction_etype_treat','WithinModel', Time, 'WithinModel', 'separatemeans');
    
% 4. Use interaction factor interaction_cat as the first variable in % multcompare, 'By' Time
tbl2=multcompare(rm2,'interaction_etype_treat','By','Time');



