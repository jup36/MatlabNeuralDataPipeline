% README for the TONIC PACKAGE

% Starting to use a new design for packages. Top-level directory will
% contain a few subfolders that organize the package. The distributable
% core of the package is a library or functions used to transform, analyze,
% and plot data.
% 
% Within this directory, the useable form of the distribution for the lab
% will contain a subfolder that holds "Analysis Workflows". This directory
% contains 2nd order functions that call upon the TONIC library to perform
% analysis by stringing together library functions into a workflow.
% 

% Algorithm for spike sorting invovles a few stages:
% 
% (0) filtering / denoising
%     acausal filtering is applied to the raw data
% (1) detection 
%     find moments at which the recording exceeds a given SNR threshold and also has a rapid recovery to a level below the SNR threshold.
% (2) extraction
%     sort the extracted events to allow largest amp events (compared to other channels) to trigger extraction of data windows from other rec channels.
% (3) alignment
%     does alignment make sense in cases of multichannel?
% (4) ts-quant
%     calculate projections onto freq or timeseries space. other scalar methods implemented as well.
% (5) sp-quant
%     calculate projections of distribution of ts-quant projections in electrode space (i.e. across channels).
% (6) id-ing
%     cluster in sp-quant or ts-quant space according to the watershed like algorithm
% 