function [L] = TNC_GM_GaussMixLogLikelihood(PI, MU, SIGMA, SIGMA_DET, SIGMA_INV, DATA, version)
    N = int32(size(DATA, 1));
    K = int32(size(DATA, 2));
%   disp(['size(DATA)=' num2str(size(DATA))]);
%   disp(['det(SIGMA)=' num2str(det(SIGMA))]);
    L_VECT = zeros(1, N);
    M = int32(size(MU, 2));  
    for i=1:N
        Xi = DATA(i,:);
        if ~isnan(Xi)
            MVG = MultiVarGaussMix(Xi, M, K, PI, MU, SIGMA, SIGMA_DET, SIGMA_INV, version);
%           if MVG < 0 || abs(imag(MVG)) > 0
%               disp(['i=' num2str(i) ' MVG=' num2str(MVG)]);
%           end
            L_VECT(i) = log(MVG);
            if isinf(L_VECT(i)) L_VECT(i)=0; end
        end
    end
    if 0
        disp('L_VECT=');
        L_VECT
        disp('SIGMA_DET=');
        SIGMA_DET
    end
    L = sum(L_VECT);

%-------------------------------------------------------------------------------

function [ Xi_contr ] = MultiVarGaussMix(Xi, M, K, PI, MU, SIGMA, SIGMA_DET, SIGMA_INV, version)
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
    %   SIGMA_INV is K x (K*M) matrix
    Xi_contr = 0;
    for m=1:M
        SIGMA1 = SIGMA(:,(K*(m-1)+1):(K*m));
        if version == 1
            Xi_contr = Xi_contr + PI(m)*mvnpdf(Xi, MU(:,m)', SIGMA1); 
        else
            Xi_contr = Xi_contr + PI(m)*TNC_GM_MultiVarGauss(Xi, MU(:,m)', ...
                                  SIGMA_INV(:,(K*(m-1)+1):(K*m)), SIGMA_DET(m));
        end
    end

