function [ X_contr ] = TNC_GM_MultiVarGauss(X, MU, SIGMA_INV, SIGMA_DET)
    K = length(X);
    X_contr = 1/sqrt((2*pi)^K * SIGMA_DET) * exp(-(X-MU)*SIGMA_INV*(X-MU)'/2);
