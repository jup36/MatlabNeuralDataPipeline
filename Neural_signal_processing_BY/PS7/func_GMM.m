function [mu, Sigma, ppi, gam, LL]=func_GMM(InitParams,Spikes)
% [mu, Sigma, ppi]=func_GMM(InitParams,Spikes)
% EM algorithm for Gaussian Mixture Model estimation
%
% xDim: data dimensionality
% zDim: number of mixture components
% N: number of data points
%
% INPUTS:
% InitParams - a 1x1 structure containing two fields
% InitParams.mu - initialization of mean vectors of GMs (xDim x zDim)
% InitParams.Sigma - initialization of covariance matrices of GMs (xDim x xDim x
% zDim)
% Spikes - input data (xDim x N)
%
% OUTPUTS:
% mu - estimated mean vectors of GMs (xDim x zDim)
% Sigma - estimated covariance matrices of GMs (xDim x xDim x zDim)
% ppi - estimated weights of GMs
% gam - estimated responsibilities of each cluster to each data point(zDim x xDim)
% LL - estimated log-likelihood at each iteration
mu = InitParams.mu; 
ppi = InitParams.pi;
K = size(mu, 2);        % number of classes
[D, N] = size(Spikes);  
Sigma = InitParams.Sigma;
const = -0.5 * D * log(2*pi);
for i = 1:100   % 100 iterations
    % === E-step ===
    logMat = nan(K, N);
    for k = 1:K
        S = Sigma(:,:,k);
        xdif = bsxfun(@minus, Spikes, mu(:,k));     % bsxfun applies element-by-element binary operation to two arrays with singleton expansion enabled
        term1 = -0.5 * sum((xdif' * inv(S)) .* xdif', 2); % N x 1
        term2 = const - 0.5 * log(det(S)) + log(ppi(k)); % scalar
        logMat(k,:) = term1' + term2;   % log(N(Xn|MuK,SigmaK)) + log(pi)
    end
    % Evaluate log P({x})
    astar = max(logMat, [], 1);
    adif = bsxfun(@minus, logMat, astar);
    nLL = log(sum(exp(adif), 1)) + astar; % 1 x N   % nLL = log likelihood?
    LL(i) = sum(nLL);
    gam = exp(bsxfun(@minus, logMat, nLL)); % K x N (responsibilities)
    gam = bsxfun(@rdivide, gam, sum(gam, 1)); % for numerical stability
    % === M-step ===
    Neff = sum(gam, 2);
    ppi = Neff' / N;
    for k = 1:K
        mu(:,k) = (Spikes * gam(k,:)') / Neff(k);
        xdif = bsxfun(@minus, Spikes, mu(:,k));
        S = bsxfun(@times, xdif, gam(k,:)) * xdif' / Neff(k);
        Sigma(:,:,k) = (S + S') / 2; % for numerical stability
    end
end
return;