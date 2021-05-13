function [rez] = leastSquareDecoder3DhTrj(s)
%This function performs 'linear least squares' decoding of 3-D hand
% trajectories from a neural spike dataset with a 'leave-a-trial-out' scheme.  
% Input: a structure 's' containing trial-by-trial neural spike data and
%   3-D hand trajectory data (position and/or velocity). 

% filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles'; 
% load(fullfile(filePath,'preprocessDecodeHTrjCtxStr_WR40_081919.mat'),'s')

if length(s.spk)==length(s.hVel)
    nTr = length(s.spk); 
else
    error('The # of trials unmatched between behavioral and neural data!')
end

%nTestTr = floor(nTr/10); 
%nTrainTr = nTr-nTestTr;

ctxId = ~s.spk_strIdx; 
strId = s.spk_strIdx; 

for k = 1:nTr
    % sort train and test trials
    
    %trialI = randperm(nTr); 
    %trainTrs = trialI(1:nTrainTr); 
    %testTrs = trialI(nTrainTr+1:end); 
    
    testTrs = false(1,length(s.spk)); 
    testTrs(k) = true;
    trainTrs = ~testTrs; 
    
    % concatanate training data
    tmpN_train = [s.spk{1,trainTrs}];
    tmpB_train = [s.hVel{1,trainTrs}];
    
    % calculate unregulated regression
    dW = pinv(tmpN_train')*tmpB_train';
    rez.dW(:,:,k) = dW; % stock decoding weights 
    
    % calculate unregulated regression separately for Ctx population
    dW_Ctx = pinv(tmpN_train(ctxId,:)')*tmpB_train';
    rez.dW_Ctx(:,:,k) = dW_Ctx; % stock decoding weights 
    
    % calculate unregulated regression separately for Str population
    dW_Str = pinv(tmpN_train(strId,:)')*tmpB_train'; 
    rez.dW_Str(:,:,k) = dW_Str; % stock decoding weights 
     
    % get estimated position for the held-out trial 
    tmpN_test = [s.spk{1,testTrs}]; % held-out neural data
    tmpB_est  = (tmpN_test'*dW)'; % estimated trj
    tmpB_shuffle_est = (tmpN_test'*dW(randperm(size(dW,1)),:))'; % estimated trj with shuffled neural-behavioral mapping
    
    tmpB_est_CtxStr = (tmpN_test(ctxId,:)'*dW(ctxId,:))'; % estimated trj by Ctx data using dW computed using both populations
    tmpB_est_StrCtx = (tmpN_test(strId,:)'*dW(strId,:))'; % estimated trj by Str data using dW computed using both populations 
    
    tmpB_est_ctx = (tmpN_test(ctxId,:)'*dW_Ctx)'; % estimated trj
    tmpB_est_str = (tmpN_test(strId,:)'*dW_Str)'; % estimated trj
   
    tmpB_test = [s.hVel{1,testTrs}]; % held-out actual hand trj data
         
    % evaluate cross-validation performance
    rez.out_cor(k) = corr(tmpB_est(:),tmpB_test(:)); 
    rez.out_corX(k) = corr(tmpB_est(1,:)',tmpB_test(1,:)'); 
    rez.out_corY(k) = corr(tmpB_est(2,:)',tmpB_test(2,:)'); 
    rez.out_corZ(k) = corr(tmpB_est(3,:)',tmpB_test(3,:)'); 
    
    rez.shuff_cor(k) = corr(tmpB_shuffle_est(:),tmpB_test(:)); 
    rez.shuff_corX(k) = corr(tmpB_shuffle_est(1,:)',tmpB_test(1,:)'); 
    rez.shuff_corY(k) = corr(tmpB_shuffle_est(2,:)',tmpB_test(2,:)'); 
    rez.shuff_corZ(k) = corr(tmpB_shuffle_est(3,:)',tmpB_test(3,:)'); 
    
    rez.out_CtxStr_cor(k) = corr(tmpB_est_CtxStr(:),tmpB_test(:)); 
    rez.out_StrCtx_cor(k) = corr(tmpB_est_StrCtx(:),tmpB_test(:)); 
   
    rez.out_ctx_cor(k) = corr(tmpB_est_ctx(:),tmpB_test(:)); 
    rez.out_str_cor(k) = corr(tmpB_est_str(:),tmpB_test(:)); 
    
    fprintf('completed fold# %d\n', k) % report fold progression
end

save(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles','leastSquareDecoder3DhVel_WR40_081919'),'rez'); % save rez

%% Assess decoded true cross validation
disp(['Mean cross-validated correlation: ' num2str(mean(rez.out_cor))]);
figure; clf;
boxplot([rez.out_cor' rez.shuff_cor' rez.out_ctx_cor' rez.out_str_cor'],'Labels',{'All cells','Shuffle','Ctx','STR'}); 
axis([0 5 -0.4 1]); ylabel('Correlation btwn decoded and actual'); box off;
print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles','Figure','boxplot_corr_decoded_actual_hVel'), '-dpdf', '-painters', '-bestfit')

%% Assess decoding across all valid trials
all_Spk = []; % ctx_Spk = []; str_Spk = []; 
rez.randp_hVel = []; 
smth_ct = 251;

dW_m = nanmean(rez.dW,3); 
dW_Ctx_m = nanmean(rez.dW_Ctx,3); 
dW_Str_m = nanmean(rez.dW_Str,3); 

% get some neural and hand velocity data
for p = randSampleToPlot% randsample(length(s.spk),10) %randperm(10)
    all_Spk = [all_Spk s.spk{1,p}];
    %ctx_Spk = [ctx_Spk s.spk{1,p}(ctxId,:)]; 
    %str_Spk = [str_Spk s.spk{1,p}(strId,:)]; 
    
    rez.randp_hVel = [rez.randp_hVel s.hVel{1,p}]; 
end

rez.all_Spk_hV = sgolayfilt(all_Spk'*dW_m, 3, 251)'; 
rez.ctx_Spk_hV = sgolayfilt(all_Spk(ctxId,:)'*dW_Ctx_m, 3, 251)'; 
rez.str_Spk_hV = sgolayfilt(all_Spk(strId,:)'*dW_Str_m, 3, 251)'; 

figure; clf; 
subplot(311); hold on; 
plot(rez.ctx_Spk_hV(1,:),'color',[0 .5 .33]); 
plot(rez.str_Spk_hV(1,:),'color',[0.33 0 .5]); 
plot(rez.randp_hVel(1,:),'k'); 
title('Actual and decoded hVel (X)')
xlabel('time')
ylabel('mm')

subplot(312); hold on; 
plot(rez.ctx_Spk_hV(2,:),'color',[0 .5 .33]); 
plot(rez.str_Spk_hV(2,:),'color',[0.33 0 .5]); 
plot(rez.randp_hVel(2,:),'k'); 
title('Actual and decoded hVel (Y)')
xlabel('time')
ylabel('mm')

subplot(313); hold on; 
plot(rez.ctx_Spk_hV(3,:),'color',[0 .5 .33]); 
plot(rez.str_Spk_hV(3,:),'color',[0.33 0 .5]); 
plot(rez.randp_hVel(3,:),'k'); 
title('Actual and decoded hVel (Z)')
xlabel('time')
ylabel('mm')
print(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles','Figure','decoded_actual_hVel_CtxStr'), '-dpdf', '-painters', '-bestfit')

%% Assess decoding trial-by-trial
nTrPlot = 1; 
figure; clf;  
for p = randsample(length(s.spk),nTrPlot)
    
    ctx_out = sgolayfilt(s.spk{1,p}(ctxId,:)'*dW_Ctx_m,3,smth_ct)';
    str_out = sgolayfilt(s.spk{1,p}(strId,:)'*dW_Str_m,3,smth_ct)';
    
    hold on; 
    plot(sgolayfilt((s.hVel{p}(1,:))',3,smth_ct)',sgolayfilt((s.hVel{p}(3,:))',3,smth_ct)','k')
    plot(ctx_out(1,:),ctx_out(3,:),'color',[0 .5 .33]);
    plot(str_out(1,:),ctx_out(3,:),'color',[0.33 0 .5]);
    hold off; 
end

end