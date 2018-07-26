function [kernel] = TNC_CreateCausalKernel(tau_rise,tau_decay,dt)

% Create an exponential, causal filter for calculating PSTHs
    kernelTMP(1,:) = zeros(1,(tau_decay.*8)./dt);
    time(1,:) = dt:dt:(tau_decay.*8);
    kernelTMP(1,:) = (1-exp(-time(1,:)./tau_rise)).*exp(-time(1,:)./tau_decay);
    kernelTMP(1,length(kernelTMP(1,:))) = 0;
    kernelTMP(1,1) = 0;

% Normalize the kernel to have unit area
    kernelTMP_INT = trapz(kernelTMP(1,:));
%     kernelTMP_INT = max(kernelTMP(1,:));
    kernel(1,:) = kernelTMP(1,:)./kernelTMP_INT;
    
    
