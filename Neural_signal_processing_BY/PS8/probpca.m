function [estParams, LL] = probpca(X, zDim, varargin)
%
% [estParams, LL] = fastfa(X, zDim, ...)
%
% Factor analysis and probabilistic PCA.
%
%   xDim: data dimensionality
%   zDim: latent dimensionality
%   N:    number of data points
%
% INPUTS:
%
% X    - data matrix (xDim x N)
% zDim - number of factors
%
% OUTPUTS:
%
% estParams.L  - factor loadings (xDim x zDim)
% estParams.Ph - diagonal of uniqueness matrix (xDim x 1)
% estParams.d  - data mean (xDim x 1)
% LL           - log likelihood at each EM iteration
%
% OPTIONAL ARGUMENTS:
%
% typ        - 'fa' (default) or 'ppca'
% tol        - stopping criterion for EM (default: 1e-8)
% cyc        - maximum number of EM iterations (default: 1e8)
% minVarFrac - fraction of overall data variance for each observed dimension
%              to set as the private variance floor.  This is used to combat
%              Heywood cases, where ML parameter learning returns one or more
%              zero private variances. (default: 0.01)
%              (See Martin & McDonald, Psychometrika, Dec 1975.)
% verbose    - logical that specifies whether to display status messages
%              (default: false)
%
% Code adapted from ffa.m by Zoubin Ghahramani.
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

  typ        = 'ppca';
  tol        = 1e-8;
  cyc        = 1e8;
  minVarFrac = 0.01;
  verbose    = false;
% assignopts(who, varargin);

  randn('state', 0);
  [xDim, N] = size(X);

  % Initialization of parameters
  cX = cov(X', 1);   % covariance of X
  if rank(cX) == xDim   % rank(A) provides an estimate of the number of linearly independent rows or columns of a matrix A
    scale = exp(2*sum(log(diag(chol(cX))))/xDim);   % geometric mean?
  else
    % cX may not be full rank because N < xDim
    fprintf('WARNING in fastfa.m: Data matrix is not full rank.\n');
    r     = rank(cX);
    e     = sort(eig(cX), 'descend');
    scale = geomean(e(1:r));
  end
  
  L     = randn(xDim,zDim)*sqrt(scale/zDim); % W, init. factor loadings (xDim x zDim)
  Ph    = diag(cX);                          % Covariance of P(X|Z), the observation noise, diagonal of uniqueness matrix (xDim x 1)
  d     = mean(X, 2);                        % data mean across N (xDim x 1)

  varFloor = minVarFrac * diag(cX); % fraction of overall data variance for each observed dimension to set as the private variance floor

  I     = eye(zDim);    % identity matrix
  const = -xDim/2*log(2*pi);    % constant term in the MVG pdf
  LLi   = 0;
  LL    = [];

  for i = 1:cyc     % maximum number of EM iterations
    % =======
    % E-step; Goal: Maximize logP(X|Params) w.r.t. Params, Where Params are L,Mu,Sigma.
    % =======
    invC = inv(L*L' + diag(Ph)); % invC = inv(Sigma), Sigma equals covariance of X (Xdim x Xdim), P(X|Params), Eq9 in the note
    beta = L' * invC;            % zDim x xDim, W*inv(Sigma)

    cX_beta = cX * beta'; % xDim x zDim, sum(x-mu)*Ez (numerator for W), (C)*W*(C^-1)
    EZZ     = I - beta * L + beta * cX_beta;    % using, E[ZZ]=Cov[Z]+E[Z]E[Z] and the distribution of P(Z|X)

    % Compute log likelihood
    LLold = LLi;    % update previous log-likelihood
    ldM   = sum(log(diag(chol(invC))));     % determinant of invC
    LLi   = N*const + N*ldM - 0.5*N*sum(sum(invC .* cX)); % current log-likelihood P(X|param)
    if verbose  % display the progress if chosen to do so
      fprintf('EM iteration %5i lik %8.1f \r', i, LLi);
    end
    LL = [LL LLi];  % update log-likelihood array

    % =======
    % M-step
    % =======
    L  = cX_beta / EZZ;     % update loading matrix(vector) using Eq12 in the note 
    Ph = diag(cX) - sum(cX_beta .* L, 2);   % Sigma (observation noise), Eq13

    if isequal(typ, 'ppca')
      Ph = mean(Ph) * ones(xDim, 1);    % ppca assumes isotropic observation noise
    end
    
    if isequal(typ, 'fa')
      % Set minimum private variance
      Ph = max(varFloor, Ph);    % FA allows non-isotropic observation noise
    end

    if i<=2
      LLbase = LLi;     % initial LLs
    elseif (LLi < LLold)
      disp('VIOLATION');    % if current LL is smaller than previous LL
    elseif ((LLi-LLbase) < (1+tol)*(LLold-LLbase))  % if the current LL only improves marginally from the previous value, BREAK the loop!
      break;
    end
  end

  if verbose
    fprintf('\n');
  end

  if any(Ph == varFloor)
    fprintf('Warning: Private variance floor used for one or more observed dimensions in FA.\n');
  end

  estParams.L  = L;
  estParams.Ph = Ph;
  estParams.d  = d;

