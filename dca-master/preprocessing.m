function preprocessing(p)
% - compute any fixed variables before optimization
% - initialize any needed quantities

num_datasets = length(p.Results.X); % added by JP to run line-by-line within the DCA.m

X = p.Results.X;  % X will change as we optimize each dim
X_orig = X;       % X_orig will remain the original X

%%% check how many dca dimensions there should be
%   for minimum number of dimensions across datasets + user input
num_dims_foreach_set = [];
for iset = 1:num_datasets
    num_dims_foreach_set = [num_dims_foreach_set size(X{iset},1)];
end
num_dca_dims = min(num_dims_foreach_set);

if (~isempty(p.Results.num_dca_dimensions))
    num_dca_dims = p.Results.num_dca_dimensions;
end

%%% compute the combined recentered matrices for given distance matrices
num_samples = size(p.Results.X{1,1},2);  % added by JP to run line-by-line within the DCA.m

R_given = zeros(num_samples); % Dimension: sample x sample 
D_given = p.Results.Ds;       % Dimension: sample x sample (distance matrix)
if (~isempty(D_given))
    for iset = 1:length(D_given) % increment the number of given dependent variables 
        H = eye(size(D_given{iset})) - 1 / size(D_given{iset},1) * ones(size(D_given{iset}));
        D_given{iset} = H * D_given{iset} * H;  % recenter D
        R_given = R_given + D_given{iset};      % combine
    end
    
    R_given = R_given / length(D_given); % the combined recentered matrices for given distance matrices (of dependent variables)
end

%%% prepare indices for column indices when subtracting off distance matrix means
if (p.Results.num_stoch_batch_samples == 0)  % full gradient descent
    col_indices = [];
    for icol = 1:num_samples
        col_indices = [col_indices icol:num_samples:num_samples^2]; % column indices
    end
else     % stochastic gradient descent
    col_indices = [];
    for icol = 1:p.Results.num_stoch_batch_samples
        col_indices = [col_indices icol:p.Results.num_stoch_batch_samples:p.Results.num_stoch_batch_samples^2];
    end
end

%%% initialize parameters
U = cell(1,num_datasets);      % cell vector, dca dimensions for each dataset
dcovs = zeros(1,num_dca_dims); % vector, dcovs for each dimension
for iset = 1:num_datasets
    U_orth{iset} = eye(size(X{iset},1));  % keeps track of the orthogonal space of u
end

%%% compute Xij (num_neurons x num_samples^2) for each dataset, where Xij = X_i - X_j
if (p.Results.num_stoch_batch_samples == 0)  % only for full grad descent
    Xij = [];
    for iset = 1:num_datasets
        % compute all combinations of differences between samples of X
        Xij{iset} = bsxfun(@minus, X{iset}, permute(X{iset}, [1 3 2])); % permute rearranges the dimensions of a matrix
        Xij{iset} = reshape(-Xij{iset}, size(X{iset},1), []); % dimension must be neuron x (sample x sample)
    end
    Xij_orig = Xij;
end

end
