function p = parse_input(X, vargs)
% parses input, and extracts name-value pairs

p = inputParser;  % creates parser object

default_D = [];   % distance matrices
default_num_iters_per_dataset = 1;
default_num_iters_foreach_dim = 30;
default_percent_increase_criterion = 0.01;  % stops when objective function increases fewer than 1% of current value
default_num_dca_dimensions = []; % number of dca dimensions to identify
default_u_0s = [];  % how to intialize the dimensions before optimization
default_num_stoch_batch_samples = 0;
default_num_samples_to_compute_stepwise_dcov = 1000;

addRequired(p, 'X'); % addRequired(p, 'Xs')
addOptional(p, 'Ds', default_D);
addParamValue(p, 'num_iters_per_dataset', default_num_iters_per_dataset);
addParamValue(p, 'num_iters_foreach_dim', default_num_iters_foreach_dim);
addParamValue(p, 'percent_increase_criterion', default_percent_increase_criterion);
addParamValue(p, 'num_dca_dimensions', default_num_dca_dimensions);
addParamValue(p, 'num_stoch_batch_samples', default_num_stoch_batch_samples);
addParamValue(p, 'num_samples_to_compute_stepwise_dcov', default_num_samples_to_compute_stepwise_dcov);
addParamValue(p, 'u_0s', default_u_0s);

% NOTE: addParamValue should be changed to addParameter...
%      addParamValue is for older matlab versions

parse(p,X,vargs{:});  % parses input to get optional parameters
% to get input, use p.Results.X, etc.
end


