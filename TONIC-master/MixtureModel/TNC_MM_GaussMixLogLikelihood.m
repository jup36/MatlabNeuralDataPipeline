%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function [L] = TNC_MM_GaussMixLogLikelihood(PI, MU, SIGMA, DATA)
    N = int32(size(DATA, 1));
    K = int32(size(DATA, 2));
%   disp(['size(DATA)=' num2str(size(DATA))]);
%   disp(['det(SIGMA)=' num2str(det(SIGMA))]);
    L_VECT = zeros(1, N);
    M = int32(size(MU, 2));  
    L = 0;
    for i=1:N
        Xi = DATA(i,:);
        if ~isnan(Xi)
            MVG = MultiVarGaussMix(Xi, M, K, PI, MU, SIGMA);
            L_VECT(i) = log(MVG);
            if isinf(L_VECT(i)) 
                L = -Inf;
                break; 
            end
        end
    end
    if ~isinf(L)
        L = sum(L_VECT);
    end

%-------------------------------------------------------------------------------

function [ Xi_contr ] = MultiVarGaussMix(Xi, M, K, PI, MU, SIGMA)
    % Compute PDF of a mixture of multivariatee normal distributions
    % for spike sorting problem
    %
    % Here,
    %   K is the number of channels/electrodes in a shank
    %   M is the number of components (gaussians) in the Gaussian mixture model
    %   Xi is a K-vector representing a data point
    %   PI is M-vector of mixing coefficients
    %   MU is K x M matrix
    %   SIGMA     is K x (K*M) matrix
    Xi_contr = 0;
    for m=1:M
        SIGMA_m = SIGMA(:,(K*(m-1)+1):(K*m));
%       SIGMA_m
%       det(SIGMA_m)
        try
            Xi_contr = Xi_contr + PI(m)*mvnpdf(Xi, MU(:,m)', SIGMA_m); 
        catch
            Xi_contr = Xi_contr;
        end
    end
