function D = get_D_uXij(u)
% compute distance matrix of (X projected to utemp)

D = squareform(pdist((u' * X)'));
end
