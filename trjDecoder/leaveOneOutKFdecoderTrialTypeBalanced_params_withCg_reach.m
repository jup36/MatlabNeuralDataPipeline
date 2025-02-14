function  [ctxRez, strRez, ctxstrRez, params] = leaveOneOutKFdecoderTrialTypeBalanced_params_withCg_reach(stateD, s, resample)
%stateD = tmpState;
%ctxD = s.dat.spkCtxR; 
%strD = s.dat.spkStrR;
%laserI = s.dat.laserIdxR; 
%resample = 20; 
%ctxI_bsc = s.ctxI; 
%strI_bsc = s.strI; 

nd = {};  % container cell for all neural data
counter = 0; 
if isfield(s.dat, 'spkCtxR')
    counter = counter + 1; 
    nd{1, counter} = 'ctx';
    nd{2, counter} = s.dat.spkCtxR;  
    nd{3, counter} = find(s.ctxI); 
    valTrI = cell2mat(cellfun(@(a) ~isempty(a), nd{2, counter}, 'un', 0)); % valid trials
    nd{4, counter} = mode(cellfun(@(a) size(a,1), nd{2, counter}(valTrI)));  
end

if isfield(s.dat, 'spkStrR')
    counter = counter + 1; 
    nd{1, counter} = 'str';
    nd{2, counter} = s.dat.spkStrR;  
    nd{3, counter} = find(s.strI); 
    nd{4, counter} = mode(cellfun(@(a) size(a,1), nd{2, counter}(valTrI)));  
end

if isfield(s.dat, 'spkCgR')
    counter = counter + 1; 
    nd{1, counter} = 'cg';
    nd{2, counter} = s.dat.spkCgR; 
    if ~exist('valTrI', 'var')
        valTrI = cell2mat(cellfun(@(a) ~isempty(a), nd{2, counter}, 'un', 0)); % valid trials
    end
    nd{4, counter} = mode(cellfun(@(a) size(a,1), nd{2, counter}(valTrI)));  
    nd{3, counter} = (1:nd{4, counter})';  
end

if isfield(s.dat, 'laserIdx')
    laserI = s.dat.laserIdxR; 
    stmTrI = cell2mat(cellfun(@(a) sum(a)>=1, laserI, 'un', 0)); % stim trials
end

trainTrN = min(sum(valTrI&~stmTrI))-1; 

for i = 1:resample %resample % repeat resampling trials
    Ns = cell2mat(nd(4, :)); % cell numbers 
    
    ctx_record = contains(cell2mat(nd(1,:)), 'ctx');
    str_record = contains(cell2mat(nd(1,:)), 'str');
    cg_record = contains(cell2mat(nd(1,:)), 'cg'); 
    
    if length(Ns) >= 3
        rs_n = min(Ns(1:2)); 
        rs_n_all = min(Ns); 
        
        
        
    elseif length(Ns) == 2 && ctx_record && str_record
        rs_n = min(Ns(1:2)); 
    elseif length(Ns) == 1 && cg_record
        rs_n = Ns; 
    
    end
    
    randCtxI = randperm(ctxN);
    randStrI = randperm(strN);
    ctxI = randCtxI(1:minN_ctx_str); % ctx cells for this iteration
    strI = randStrI(1:minN_ctx_str); % str cells for this iteration
    
    for r = 1:size(ctxD,1) % row: trials
        for c = 1:size(ctxD,2) % column: position/torque pairs
            trainI = true(size(valTrI,1),size(valTrI,2));
            if valTrI(r,c)
                trainI(r,c) = false; % to leave one trial out as a test trial
                testI = ~trainI;     % index for the one test trial left out
                valTrainI = trainI & valTrI & ~stmTrI; % valid train trial index (exclude stim trials)
                
                % get current train trials by resampling (to include the same # of trials for each trial type)
                trainTrs = cell(1,size(trainI,2));
                trainState = []; % current training data states (hand position: row 1-3, velocity: row 4-6)
                trainCtx = [];
                trainStr = [];
                trainCtxStr = []; 
                
                for cc = 1:size(valTrainI,2) % sample trials from each trial type
                    tempValTr = find(valTrainI(:,cc));
                    tempValTrRand = tempValTr(randperm(length(tempValTr)));
                    trainTrs{1,cc} =tempValTrRand(1:trainTrN); %tempValTrRand(1:trainTrN); % take the set # of randomized trials from each trial type (column)
                    trainState = [trainState; stateD(trainTrs{1,cc},cc)]; % concatanate randomly selected trials from each type to construct the train data matrix
                    trainCtx = [trainCtx; cellfun(@(a) a(ctxI,:), ctxD(trainTrs{1,cc},cc), 'un', 0)]; % ctx spike mat with cells resampled to match # of cells
                    trainStr = [trainStr; cellfun(@(a) a(strI,:), strD(trainTrs{1,cc},cc), 'un', 0)]; % str spike mat with cells resampled to match # of cells
                    trainCtxStr = [trainCtxStr; cellfun(@(a,b) [a(ctxI,:); b(strI,:)], ctxD(trainTrs{1,cc},cc), strD(trainTrs{1,cc},cc), 'un', 0)]; 
                end
                clearvars cc
                [NTrain,NTrType] = size(trainState);
                
                %% %%%%%%%%%%% Training phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                interStateM = 0;    % inter-state matrix
                intraStateM = 0;    % intra-state matrix
                obsStateCtxM = 0;   % observed Ctx state matrix
                obsStateStrM = 0;   % observed Str state matrix
                obsStateCtxStrM = 0;% observed CtxStr state matrix 
                obsIntraStateM = 0; % observed intra-state matrix
                count = 0;
                countPool = 0;
                for t = 1:size(trainState,1) % # of trial
                    % for parameter A, state model describes how state evolves over time
                    stateZ1 = trainState{t}(:,2:end); % Zt
                    stateZ2 = trainState{t}(:,1:end-1); % Zt-1
                    tmp1 = stateZ1*stateZ2'; % Zt*Zt-1' (6-by-6 matrix)
                    interStateM = interStateM+tmp1; % sum Zt*Zt-1'
                    tmp2 = stateZ2*stateZ2'; % Zt-1*Zt-1' (6-by-6 matrix)
                    intraStateM = intraStateM+tmp2; % sum Zt-1*Zt-1'
                    % for parameter C, observation model describes how observation relates to the state
                    stateZ = trainState{t}; % state (position, velocity)
                    obsCtx = trainCtx{t}; % ctx spike data
                    obsStr = trainStr{t}; % str spike data
                    obsCtxStr = trainCtxStr{t}; % ctx and str spike data 
                    tmp3 = obsCtx*stateZ'; % Xt*Zt' (Ctx cell-by-state mat)
                    obsStateCtxM = obsStateCtxM+tmp3; % sum Ctx observation-state matrix
                    tmp4 = obsStr*stateZ'; % Xt*Zt' (Str cell-by-state mat)
                    obsStateStrM = obsStateStrM+tmp4; % sum Str observation-state matrix
                    tmp9 = obsCtxStr*stateZ'; % Xt*Zt' (CtxStr cell-by-state mat)
                    obsStateCtxStrM = obsStateCtxStrM+tmp9; % sum CtxStr observation-state matrix
                    tmp5 = stateZ*stateZ'; % Zt*Zt'
                    obsIntraStateM = obsIntraStateM+tmp5; % sum Zt*Zt'
                    % for parameter Pi and V
                    count=count+1;
                    poolZStart(:,count) = trainState{t}(:,1);
                end
                clearvars t tt
                A=(interStateM/count)/(intraStateM/count); % Slope for the state model, which defines how state evolves over time
                C_ctx=(obsStateCtxM/count)/(obsIntraStateM/count); % Slope for the observation model for ctx spikes, which defines how observation relates to the state
                C_str=(obsStateStrM/count)/(obsIntraStateM/count); % Slope for the observation model for str spikes
                C_ctxstr = (obsStateCtxStrM/count)/(obsIntraStateM/count); % Slope for the observation model for ctx str spikes 
                Pi=mean(poolZStart,2); % Mean for the initial state (sample mean)
                V=cov(poolZStart');    % Covariance for the initial state (sample covariance)
                clearvars poolZStart
                
                % Fit Q and R (variance)
                sumQM = 0;
                sumRM_ctx = 0;
                sumRM_str = 0;
                sumRM_ctxstr = 0; 
                count = 0;
                for t = 1:size(trainState,1) % # of trial
                    % for parameter Q, covariance of the state model
                    stateZ1 = trainState{t}(:,2:end); % Zt
                    stateZ2 = trainState{t}(:,1:end-1); % Zt-1
                    tmp6 = stateZ1-A*stateZ2; % Zt-A*Zt-1
                    tmp6 = tmp6*tmp6'; % (Zt-A*Zt-1)(Zt-A*Zt-1)'
                    count = count + size(stateZ1,2);
                    sumQM = sumQM + tmp6; % Sum (Zt-A*Zt-1)(Zt-A*Zt-1)'
                    % for parameter R of Ctx spikes, covariance of the observation model
                    stateZ = trainState{t}; % Zt
                    obsCtx = trainCtx{t}; % ctx spike data
                    tmp7 = obsCtx-C_ctx*stateZ; % Xt_ctx-C_ctx*Zt
                    tmp7 = tmp7*tmp7'; % (Xt_ctx-C_ctx*Zt)(Xt_ctx-C_ctx*Zt)
                    sumRM_ctx = sumRM_ctx+tmp7; % Sum (Xt_ctx-C_ctx*Zt)(Xt_ctx-C_ctx*Zt)
                    % for parameter R of Str spikes, covariance of the observation model
                    obsStr = trainStr{t}; % str spike data
                    tmp8 = obsStr-C_str*stateZ; % Xt_str-C_str*Zt
                    tmp8 = tmp8*tmp8'; % (Xt_str-C_str*Zt)(Xt_str-C_str*Zt)
                    sumRM_str = sumRM_str+tmp8; % Sum (Xt_str-C_str*Zt)(Xt_str-C_str*Zt)
                    % for parameter R of CtxStr spikes, covariance of the observation model
                    obsCtxStr = trainCtxStr{t}; % ctx str spike data
                    tmp10 = obsCtxStr-C_ctxstr*stateZ; % Xt_ctxstr-C_ctxstr*Zt
                    tmp10 = tmp10*tmp10'; % (Xt_ctxstr-C_ctxstr*Zt)(Xt_ctxstr-C_ctxstr*Zt)
                    sumRM_ctxstr = sumRM_ctxstr+tmp10; % Sum (Xt_ctxstr-C_ctxstr*Zt)(Xt_ctxstr-C_ctxstr*Zt)    
                end
                clearvars t tt
                Q = sumQM/count; % count 1/(t-1), count should correspond to 2 to T
                R_ctx = sumRM_ctx/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
                R_str = sumRM_str/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
                R_ctxstr = sumRM_ctxstr/(count+NTrain*NTrType); 
                %% %%%%%%%%%%%%%% Test phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                curTrLength = size(stateD{testI},2); % position, velocity in 20ms bins
                % Initialization
                mu_ctx=Pi;   % sample mean
                sigma_ctx=V; % sample covariance
                mu_str=Pi;   % sample mean
                sigma_str=V; % sample covariance
                mu_ctxstr = Pi; % sample mean
                sigma_ctxstr = V; % sample covariance
                
                clearvars estStateMean_*;
                clearvars estStateCov_*;
                valCellICtx = ~(sum(C_ctx,2)==0)&~isnan(sum(C_ctx,2));
                valCellIStr = ~(sum(C_str,2)==0)&~isnan(sum(C_str,2));
                valCellICtxStr = ~(sum(C_ctxstr,2)==0)&~isnan(sum(C_ctxstr,2)); %~(sum(C_ctxstr,2)==0)&~isnan(sum(C_ctxstr,2));
                
                C_ctxVal = C_ctx(valCellICtx,:); % to drop the cells with no spike at all from the mapping C, if not it leads to singular matrix warning
                R_ctxVal = R_ctx(valCellICtx,valCellICtx); 
                C_strVal = C_str(valCellIStr,:); 
                R_strVal = R_str(valCellIStr,valCellIStr); 
                C_ctxstrVal = C_ctxstr(valCellICtxStr,:);
                R_ctxstrVal = R_ctxstr(valCellICtxStr,valCellICtxStr);
                
                tempC_ctx = nan(length(ctxI_bsc),size(C_ctxVal,2)); 
                tempC_ctx(ctxI_bscN(ctxI(valCellICtx)),:) = C_ctxVal; 
                params.C_ctxVal{r,c,i} = tempC_ctx; 
                
                tempC_str = nan(length(strI_bsc),size(C_strVal,2)); 
                tempC_str(strI_bscN(strI(valCellIStr)),:) = C_strVal; 
                params.C_strVal{r,c,i} = tempC_str; 
                
                tempC_CtxStr = nan(length(ctxI_bsc),size(C_ctxstrVal,2)); 
                tempC_CtxStrI = [ctxI_bscN(ctxI); strI_bscN(strI)]; 
                tempC_CtxStrI(valCellICtxStr); 
                tempC_CtxStr(tempC_CtxStrI(valCellICtxStr),:)=C_ctxstrVal; 
                params.C_ctxstrVal{r,c,i} = tempC_CtxStr; 
                                
                for b = 1:curTrLength % # of timebins
                    % One-step prediction by ctx and str data separately
                    mu_ctx = A*mu_ctx; % A is the coefficient matrix for the state model that maps the previous states to current states
                    sigma_ctx = A*sigma_ctx*A'+Q; % Q is the covariance of the state model
                    mu_str = A*mu_str;
                    sigma_str = A*sigma_str*A'+Q;
                    mu_ctxstr = A*mu_ctxstr; 
                    sigma_ctxstr = A*sigma_ctxstr*A'+Q;
                    
                    % compute the Kalman gain (needs to separately computed for ctx and str)
                    K_ctx = sigma_ctx*C_ctxVal'/(C_ctxVal*sigma_ctx*C_ctxVal'+R_ctxVal); % *inv(C_ctx*sigma*C_ctx'+R_ctx); % Kalman gain for ctx observation model
                    K_str = sigma_str*C_strVal'/(C_strVal*sigma_str*C_strVal'+R_strVal); % *inv(C_str*sigma*C_str'+R_str); % Kalman gain for str observation model
                    K_ctxstr = sigma_ctxstr*C_ctxstrVal'/(C_ctxstrVal*sigma_ctxstr*C_ctxstrVal'+R_ctxstrVal); 
                    
                    % update the state by ctx observation
                    curObsX_ctx = ctxD{testI}(ctxI,b); % take the current ctx spikes bin-by-bin
                    curObsX_ctx = curObsX_ctx(valCellICtx,1); 
                    
                    mu_ctx = mu_ctx + K_ctx*(curObsX_ctx-C_ctxVal*mu_ctx);
                    sigma_ctx = sigma_ctx - K_ctx*C_ctxVal*sigma_ctx;
                    estStateMean_ctx(:,b)=mu_ctx;
                    estStateCov_ctx(:,:,b)=sigma_ctx;
                    
                    % update the state by str observation
                    curObsX_str = strD{testI}(strI,b); % take the current str spikes bin-by-bin
                    curObsX_str = curObsX_str(valCellIStr,1); 
                    
                    mu_str = mu_str + K_str*(curObsX_str-C_strVal*mu_str);
                    sigma_str = sigma_str - K_str*C_strVal*sigma_str;
                    estStateMean_str(:,b)=mu_str;
                    estStateCov_str(:,:,b)=sigma_str;
                    
                    % update the state by ctxstr observation
                    curObsX_ctxstr = [ctxD{testI}(ctxI,b);strD{testI}(strI,b)]; % take the current ctx spikes bin-by-bin
                    curObsX_ctxstr = curObsX_ctxstr(valCellICtxStr,1); %curObsX_ctxstr([valCellICtx;valCellIStr],1); 
                    
                    mu_ctxstr = mu_ctxstr + K_ctxstr*(curObsX_ctxstr-C_ctxstrVal*mu_ctxstr);
                    sigma_ctxstr = sigma_ctxstr - K_ctxstr*C_ctxstrVal*sigma_ctxstr;
                    estStateMean_ctxstr(:,b)=mu_ctxstr;
                    estStateCov_ctxstr(:,:,b)=sigma_ctxstr;             
                end
                clearvars b
                estStateCtxMean{r,c,i} = estStateMean_ctx;
                estStateCtxCov{r,c,i} = estStateCov_ctx;
                
                estStateStrMean{r,c,i} = estStateMean_str;
                estStateStrCov{r,c,i} = estStateCov_str;
                
                estStateCtxStrMean{r,c,i} = estStateMean_ctxstr;
                estStateCtxStrCov{r,c,i} = estStateCov_ctxstr;               
                %clearvars estStateMean_* estStateCov_*              
            end
        end
    end
    fprintf('finished iteration# %d\n', i)
end
clearvars r c

[ctxRez.est,ctxRez.estCut,ctxRez.estInt] = getRepresentativeStateEst(stateD,estStateCtxMean,valTrI);
[strRez.est,strRez.estCut,strRez.estInt] = getRepresentativeStateEst(stateD,estStateStrMean,valTrI);
[ctxstrRez.est,ctxstrRez.estCut,ctxstrRez.estInt] = getRepresentativeStateEst(stateD,estStateCtxStrMean,valTrI);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    function [representState,representStateCut,representIntState] = getRepresentativeStateEst(refState,estState,valTrId)
        %refState = stateD;
        %estState = estStateCtxMean;
        %valTrId = valTrI;
        nKv = mode(cell2mat(cellfun(@(a) size(a,1),refState(valTrId),'un',0)));
        for rr = 1:size(refState,1)
            for cl = 1:size(refState,2)
                if valTrId(rr,cl)
                    distToRef = squeeze(cell2mat(cellfun(@(a) sum(sum(sqrt((a-refState{rr,cl}).^2))),estState(rr,cl,:),'un',0))); % average across ctx resampled trials
                    representState{rr,cl} = estState{rr,cl,find(distToRef==min(distToRef),1,'first')};
                    tmpState = representState{rr,cl};
                    tmpCutState = nan(nKv,50);
                    if size(tmpState,2)>=5
                        tmpCutState(:,1:min(size(tmpState,2),50)) = tmpState(:,1:min(size(tmpState,2),50));
                        representStateCut{rr,cl} = tmpCutState;
                        representIntState{rr,cl} = intm(tmpState,50); % original estimated trajectory with interpolation                   
                    end
                end
            end
        end
        clearvars r c
    end

end
