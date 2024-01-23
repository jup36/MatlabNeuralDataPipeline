%% load dPrm data 
load(fullfile('/Volumes/Extreme SSD/js2p0/collectData','dPrime_CtxStr_vec_norm_collectRez'),'dPrmRez')
isStrMat = cell2mat({dPrmRez(:).isStr})';
dPrmRezCtx = dPrmRez(~isStrMat);
dPrmRezStr = dPrmRez(isStrMat);
clearvars dPrmRez

load(fullfile('/Volumes/Extreme SSD/js2p0/collectData','dPrime_Cg_vec_norm_collectRez'),'dPrmRez')
dPrmRezCg = dPrmRez; clearvars dPrmRez 

dPrmRez = [dPrmRezCtx, dPrmRezStr, dPrmRezCg]; 

% calculate the population size (# of units) normalization factor 
total = length(dPrmRezCtx) + length(dPrmRezStr) + length(dPrmRezCg); % the total number of neurons across all three areas
total_over_str = total/length(dPrmRezStr); % striatum has the largest number of units
norm_factor_str = total/length(dPrmRezStr)/total_over_str; % normalization factor for striatum  
norm_factor_ctx = total/length(dPrmRezCtx)/total_over_str; % normalization factor for cortex
norm_factor_cg = total/length(dPrmRezCg)/total_over_str; % normalization factor for cg

[vf_ctx, vf_ctx_weighted] = vector_fields(dPrmRezCtx, norm_factor_ctx); 
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/dPrm', 'dPrm_vector_fields_ctx'), '-dpdf', '-vector')

[vf_str, vf_str_weighted] = vector_fields(dPrmRezStr, norm_factor_str); 
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/dPrm', 'dPrm_vector_fields_str'), '-dpdf', '-vector')

[vf_cg, vf_cg_weighted] = vector_fields(dPrmRezCg, norm_factor_cg); 
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/dPrm', 'dPrm_vector_fields_Cg'), '-dpdf', '-vector')

[vf, vf_weighted] = vector_fields(dPrmRez, 1); 
%print(fullfile('/Volumes/Extreme SSD/js2p0/collectData/collectFigure/dPrm', 'dPrm_vector_fields'), '-dpdf', '-vector')

