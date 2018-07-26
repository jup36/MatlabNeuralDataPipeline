function check_input(p)
% make sure user inputs correct formats

%%% check X
if (~iscell(p.Results.X) || size(p.Results.X,1) > 1 && size(p.Results.X,2) > 1) % check if X is a cell vector
    error('Xs (1 x num_datasets) should be a cell array, where Xs{iset} is (num_variables x num_datapoints)');
end

num_datasets = length(p.Results.X);

[num_vars, num_samples] = cellfun(@size, p.Results.X);
if (length(unique(num_samples)) ~= 1)  % there should only be one sample size
    error('Dataset(s) in Xs do not contain the same number of samples. Xs{iset} (num_variables x num_samples), where num_samples is the same for each dataset (but num_variables can be different).');
end
num_samples = size(p.Results.X{1},2);

isnan_found = false;
for iset = 1:num_datasets
    isnan_found = isnan_found | any(any(isnan(p.Results.X{iset})));
end
if (isnan_found == true)
    error('Dataset(s) in Xs contain NaNs. Remove samples with NaNs.');
end

%%% check Ds
if (~isempty(p.Results.Ds) && (~iscell(p.Results.Ds) || size(p.Results.Ds, 1) > 1 && size(p.Results.Ds, 2) > 1))
    error('Ds should either be empty (Ds = []) or a cell vector');
end

if (length(p.Results.Ds) > 0)
    [num_samples1, num_samples2] = cellfun(@size, p.Results.Ds);
    if (length(unique([num_samples1, num_samples2])) ~= 1)
        error('Dataset(s) in Ds do not contain the same number of samples. Ds{iset} (num_samples x num_samples) are distance matrices, where num_samples is the same for each dataset.');
    end
    
    isnan_found = false;
    for iset = 1:length(p.Results.Ds)
        isnan_found = isnan_found | any(any(isnan(p.Results.Ds{iset})));
    end
    if (isnan_found == true)
        error('Dataset(s) in Ds contain NaNs. Remove samples with NaNs.');
    end
    
    isneg_found = false;
    for iset = 1:length(p.Results.Ds)
        isneg_found = isneg_found | any(any(p.Results.Ds{iset} < 0));
    end
    if (isneg_found == true)
        error('Dataset(s) in Ds contain negative values. Ds{iset} is a distance matrix with nonnegative values.');
    end
end

%%% check that X and D have more than just one dataset combined
if (length(p.Results.X) + length(p.Results.Ds) <= 1)
    error('Not enough datasets in Xs and Ds.  The number of datasets (including given distance matrices) should be at least two.');
end

%%% check that X and D have the same number of samples
if (~isempty(p.Results.Ds) && num_samples ~= num_samples1)
    error('Dataset(s) in Xs do not have the same number of samples as those in Ds.  They should be the same.');
end

%%% check num_dca_dimensions
if (p.Results.num_dca_dimensions > min(num_vars))
    error(sprintf('"num_dca_dimensions" must be less than or equal to %d, the minimum number of variables across datasets.', min(num_vars)));
end
end
