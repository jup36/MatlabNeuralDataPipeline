dPrmRez_ctx_str = load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'dPrime_CtxStr_vec_norm_collectRez.mat'), 'dPrmRez'); 
dPrmRez_ctx_str = dPrmRez_ctx_str.('dPrmRez'); 
clearvars dPrmRez

isStrI = [dPrmRez_ctx_str.isStr]; 
sigI = [dPrmRez_ctx_str.sigI]; 

ctx_sig_frac = sum(~isStrI & sigI)/sum(~isStrI); 
str_sig_frac = sum(isStrI & sigI)/sum(isStrI); 


dPrmRez_cg = load(fullfile('/Volumes/Extreme SSD/js2p0/collectData', 'dPrime_Cg_vec_norm_collectRez.mat'), 'dPrmRez'); 
dPrmRez_cg = dPrmRez_cg.('dPrmRez'); 
clearvars dPrmRez

sigI_cg = [dPrmRez_cg.sigI]; 
cg_sig_frac = sum(sigI_cg)/length(sigI_cg); 