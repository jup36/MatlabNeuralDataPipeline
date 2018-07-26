function R = get_recentered_matrix(u, X)
    % computes the recentered distance matrix for each dataset
    % u: (num_variables x 1),  weight vector
    % X: (num_variables x num_datapoints), one dataset

    % compute distance matrix of X projected onto u
        D = squareform(pdist((u' * X)'));
        
    % now recenter it
        H = eye(size(D)) - 1/size(D,1) * ones(size(D));
        R = H * D * H;  % recenters distance matrix
end