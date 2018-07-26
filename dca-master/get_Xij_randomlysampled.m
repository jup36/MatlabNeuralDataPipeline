function Xij_sampled = get_Xij_randomlysampled(X)
% computes Xij for a subsample   (for stochastic gradient descent)

    % compute all combinations of differences between samples of X
    Xij_sampled = bsxfun(@minus, X, permute(X, [1 3 2]));
    Xij_sampled = reshape(-Xij_sampled, size(X,1), []);
end
