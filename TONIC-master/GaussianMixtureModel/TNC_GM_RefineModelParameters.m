%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function [MU0,SIGMA0,EPS,stop] = ...
         TNC_GM_RefineModelParameters(MU, SIGMA, options, stop)
    K = size(SIGMA, 1);
    M = size(MU, 2);

    SIGMA_REF = zeros(K, M*K);
    SIGMA_INV = zeros(K, M*K);
    SIGMA_DET = zeros(1, M);
    ALPHA     = zeros(K, M);
    MU0       = zeros(1, M);
    SIGMA0    = zeros(1, M);
    EPS       = zeros(1, K);

    % Determine MU0 and ALPHA
    for m=1:M
        MU0(m) = mean(MU(:,m));
        denom = MU0(m);
        if denom == 0
            denom = 1;
        end
        ALPHA(:, m) = MU(:,m)/denom;  
    end

    if options.verbose
        ALPHA 
    end

    % Determine SIGMA0
    for m=1:M

        % Determine SIGMA0
        N = 0;
        SUM_SIGMA = 0;
        SUM_ALPHA = 0;
        SIGMA1 = SIGMA(:,(K*(m-1)+1):(K*m));
        for i=1:K
            for j = 1:K
                if i ~= j
                    if SIGMA1(i, j) > 0
                        SUM_SIGMA = SUM_SIGMA + SIGMA1(i, j);
                        SUM_ALPHA = SUM_ALPHA + ALPHA(i, m)*ALPHA(j, m);
                        N = N+1;
                    end
                end
            end
        end
        if SUM_ALPHA == 0
            SIGMA0(m) = SUM_SIGMA;
        else
            SIGMA0(m) = SUM_SIGMA/SUM_ALPHA;
        end
    end

    % Compute EPS
    for k=1:K
        EPS(k) = 0;
        for m=1:M
            SIGMA_term = SIGMA(k,K*(m-1)+k);
            ALPHA_term = SIGMA0(m)*ALPHA(k,m)*ALPHA(k,m);
            if  EPS(k) < SIGMA_term - ALPHA_term
                EPS(k) = SIGMA_term - ALPHA_term;
            end
        end
%       if EPS(k) == 0
%           EPS(k) = options.min_eps;
%       end
    end
    if options.verbose
        disp('In TNC_GM_RefineModelParameters.m:');
        EPS
    end
    
    % Detrmine SIGMA_REF, SIGMA_DET and SIGMA_INV
    if options.model == 2
        for m=1:M
            SIGMA_REF(:,(K*(m-1)+1):(K*m)) = ...
                get_refined_covariance_matrix(K, SIGMA0(m), ALPHA(:,m)', EPS);
            SIGMA_DET(m) = get_refined_determinant(SIGMA0(m), ALPHA(:,m)', EPS);
            SIGMA_INV(:,(K*(m-1)+1):(K*m)) = inv(SIGMA_REF(:,(K*(m-1)+1):(K*m)));
        end
    else
        SIGMA_REF = SIGMA;
        for m=1:M
            SIGMA1 = SIGMA(:,(K*(m-1)+1):(K*m));
            SIGMA_DET(m) = det(SIGMA1);
            % Regularize degenerate matrix
            if SIGMA_DET(m) <= 0    
                break;
                stop = 1;
                disp(' ');
                disp(['NOTE: SIGMA with DET(' num2str(m) ')=' num2str(SIGMA_DET(m))]);
                return;
            end
            SIGMA_INV(:,(K*(m-1)+1):(K*m)) = inv(SIGMA1);
        end
    end

%-------------------------------------------------------------------------------

function [ SIGMA_MATR ] = get_refined_covariance_matrix(K, SIGMA0, ALPHA1, EPS)
    SIGMA_MATR = zeros(K, K);
    for k=1:K
        for j=1:K
            SIGMA_MATR(k,j) = SIGMA0*ALPHA1(k)*ALPHA1(j);
            if k == j
                SIGMA_MATR(k,j) = SIGMA_MATR(k,j) + max(EPS(k), 1.e-9);
            end
        end
    end

% -------------------------------------------------------------------------------

function SIGMA_DET = get_refined_determinant(SIGMA0, ALPHA1, EPS)
    K = size(EPS, 2);
    term1 = 1;
    for k=1:K
        term1 = term1 * EPS(k);
    end
    term2 = 0;
    for k=1:K
        prod2 = ALPHA1(k)^2;
        for i=1:K
            if i ~= k
                prod2 = prod2 * EPS(i);
            end
        end
        term2 = term2 + prod2;
    end
    SIGMA_DET = term1 + term2;

