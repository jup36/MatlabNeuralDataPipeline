function [mu, Sigma, ppi, gam, LL]=func_GMM_crossval(tempInitParams, train, test)
% [mu, Sigma, ppi]=func_GMM(tempInitParams,Spikes)
% EM algorithm for Gaussian Mixture Model estimation
%
% xDim: data dimensionality
% zDim: number of mixture components
% N: number of data points
%
% INPUTS:
% tempInitParams - a 1x1 structure containing two fields
% tempInitParams.mu - initialization of mean vectors of GMs (xDim x zDim)
% tempInitParams.Sigma - initialization of covariance matrices of GMs (xDim x xDim x
% zDim)
% train - train data: xDim x n of (J-1) folds x folds
% test - test data: xDim x n of 1 fold x folds
%
% OUTPUTS:
% mu - estimated mean vectors of GMs (xDim x zDim)
% Sigma - estimated covariance matrices of GMs (xDim x xDim x zDim)
% ppi - estimated weights of GMs
% gam - estimated responsibilities of each cluster to each data point(zDim x xDim)
% LL - estimated log-likelihood at each iteration

fold = size(train,1);   % cross-validation folds equals to the number of train or test data sets 
mu = tempInitParams.mu; 
ppi = tempInitParams.pi;
K = size(mu, 2);        % number of clusters
[D, N] = size(train{1,1});      % get dimension D, the number of trials in each train set N 
Sigma = tempInitParams.Sigma;
const = -0.5 * D * log(2*pi);

for f = 1:fold      % the # of cross-validation 
    
    % Train phase
    trD = train{f,1};   % train data 
    for i = 1:100   % 100 iterations
        % === E-step ===
        logMat = nan(K, N);
        for k = 1:K     % the # of classes (clusters)
            S = Sigma(:,:,k);
            xdif = bsxfun(@minus, trD, mu(:,k));     % bsxfun applies element-by-element binary operation to two arrays with singleton expansion enabled
            term1 = -0.5 * sum((xdif' * inv(S)) .* xdif', 2); % N x 1
            term2 = const - 0.5 * log(det(S)) + log(ppi(k));  % scalar
            logMat(k,:) = term1' + term2;   % log(N(Xn|MuK,SigmaK)) + log(pi)
        end
        % Evaluate log P({x})
        astar = max(logMat, [], 1);
        adif = bsxfun(@minus, logMat, astar);
        trnLL = log(sum(exp(adif), 1)) + astar; % 1 x N, nLL = log likelihood, sum across the clusters
        trLL(i) = sum(trnLL);   % sum across the trials
        trgam = exp(bsxfun(@minus, logMat, trnLL)); % K x N (responsibilities)
        trgam = bsxfun(@rdivide, trgam, sum(trgam, 1)); % for numerical stability
        
        % === M-step ===
        Neff = sum(trgam, 2);
        ppi = Neff' / N;
        for k = 1:K
            mu(:,k) = (trD * trgam(k,:)') / Neff(k);
            xdif = bsxfun(@minus, trD, mu(:,k));
            S = bsxfun(@times, xdif, trgam(k,:)) * xdif' / Neff(k);
            Sigma(:,:,k) = (S + S') / 2; % for numerical stability
        end
    end
    clearvars i astar adif nLL
    
    % Test phase
    teD = test{f,1};    % test data
    telogMat = nan(K, size(test{1,1},2));       % cluster by data points
    for k = 1:K     % the # of classes (cluesters)
        S = Sigma(:,:,k);
        xdif = bsxfun(@minus, teD, mu(:,k));     % bsxfun applies element-by-element binary operation to two arrays with singleton expansion enabled
        term1 = -0.5 * sum((xdif' * inv(S)) .* xdif', 2); % N x 1
        term2 = const - 0.5 * log(det(S)) + log(ppi(k)); % scalar
        telogMat(k,:) = term1' + term2;   % log(N(Xn|MuK,SigmaK)) + log(pi)
    end
    
    % Evaluate log P({x}) for the test data
    astar = max(telogMat, [], 1);
    adif = bsxfun(@minus, telogMat, astar);
    tenLL = log(sum(exp(adif), 1)) + astar; % 1 x N   nLL = log likelihood, sum across the clusters
    teLL(f) = sum(tenLL);     % log likelihood, sum across the trials
    tegam = exp(bsxfun(@minus, telogMat, tenLL)); % K x N (responsibilities)
    tegam = bsxfun(@rdivide, tegam, sum(tegam, 1)); % for numerical stability
end

% The CV likelihoods for K=NumClusterList(clusterIX)
totLL = sum(teLL);     % sum the log likelihoods across the cross-validation folds

return;