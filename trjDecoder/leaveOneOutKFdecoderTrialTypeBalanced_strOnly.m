function  [strRez] = leaveOneOutKFdecoderTrialTypeBalanced_strOnly(stateD,strD,resample)
%stateD = s.dat.stateR;
%ctxD = s.dat.spkCtxR; 
%strD = s.dat.spkStrR;
%laserI = s.dat.laserIdxR; 
%resample = 20; 

valTrI0 = cell2mat(cellfun(@(a) ~isempty(a), strD, 'un', 0)); % valid trials
valTrI1 = cell2mat(cellfun(@(a) sum(sum(~isnan(a)))>0, stateD, 'un', 0));
valTrI = valTrI0 & valTrI1; 

%stmTrI = cell2mat(cellfun(@(a) sum(a)>=1, laserI, 'un', 0)); % stim trials
strN = mode(cellfun(@(a) size(a,1), strD(valTrI)));   
minN = min(strN); 
trainTrN = min(sum(valTrI))-1; 

for i = 1:resample %resample % repeat resampling trials
    randStrI = randperm(strN);
    strI = randStrI(1:minN); % str cells for this iteration
    
    for r = 1:size(strD,1) % row: trials
        for c = 1:size(strD,2) % column: position/torque pairs
            trainI = true(size(valTrI,1),size(valTrI,2));
            if valTrI(r,c)
                trainI(r,c) = false; % to leave one trial out as a test trial
                testI = ~trainI;     % index for the one test trial left out
                valTrainI = trainI & valTrI; % valid train trial index (exclude stim trials)
                
                % get current train trials by resampling (to include the same # of trials for each trial type)
                trainTrs = cell(1,size(trainI,2));
                trainState = []; % current training data states (hand position: row 1-3, velocity: row 4-6)
                trainStr = [];
                
                for cc = 1:size(valTrainI,2) % sample trials from each trial type
                    tempValTr = find(valTrainI(:,cc));
                    tempValTrRand = tempValTr(randperm(length(tempValTr)));
                    trainTrs{1,cc} =tempValTrRand(1:trainTrN); %tempValTrRand(1:trainTrN); % take the set # of randomized trials from each trial type (column)
                    trainState = [trainState; stateD(trainTrs{1,cc},cc)]; % concatanate randomly selected trials from each type to construct the train data matrix
                    trainStr = [trainStr; cellfun(@(a) a(strI,:), strD(trainTrs{1,cc},cc), 'un', 0)]; % str spike mat with cells resampled to match # of cells
                end
                clearvars cc
                [NTrain,NTrType] = size(trainState);
                
                %% %%%%%%%%%%% Training phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                interStateM = 0;    % inter-state matrix
                intraStateM = 0;    % intra-state matrix
                obsStateStrM = 0;   % observed Str state matrix
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
                    obsStr = trainStr{t}; % str spike data
                    tmp4 = obsStr*stateZ'; % Xt*Zt' (Str cell-by-state mat)
                    obsStateStrM = obsStateStrM+tmp4; % sum Str observation-state matrix
                    tmp5 = stateZ*stateZ'; % Zt*Zt'
                    obsIntraStateM = obsIntraStateM+tmp5; % sum Zt*Zt'
                    % for parameter Pi and V
                    count=count+1;
                    poolZStart(:,count) = trainState{t}(:,1);
                end
                clearvars t tt
                A=(interStateM/count)/(intraStateM/count); % Slope for the state model, which defines how state evolves over time
                C_str=(obsStateStrM/count)/(obsIntraStateM/count); % Slope for the observation model for str spikes
                Pi=mean(poolZStart,2); % Mean for the initial state (sample mean)
                V=cov(poolZStart');    % Covariance for the initial state (sample covariance)
                clearvars poolZStart
                
                % Fit Q and R (variance)
                sumQM = 0;
                sumRM_str = 0;
                count = 0;
                for t = 1:size(trainState,1) % # of trial
                    % for parameter Q, covariance of the state model
                    stateZ1 = trainState{t}(:,2:end); % Zt
                    stateZ2 = trainState{t}(:,1:end-1); % Zt-1
                    tmp6 = stateZ1-A*stateZ2; % Zt-A*Zt-1
                    tmp6 = tmp6*tmp6'; % (Zt-A*Zt-1)(Zt-A*Zt-1)'
                    count = count + size(stateZ1,2);
                    sumQM = sumQM + tmp6; % Sum (Zt-A*Zt-1)(Zt-A*Zt-1)'
                    stateZ = trainState{t}; % Zt
                    % for parameter R of Str spikes, covariance of the observation model
                    obsStr = trainStr{t}; % str spike data
                    tmp8 = obsStr-C_str*stateZ; % Xt_str-C_str*Zt
                    tmp8 = tmp8*tmp8'; % (Xt_str-C_str*Zt)(Xt_str-C_str*Zt)
                    sumRM_str = sumRM_str+tmp8; % Sum (Xt_str-C_str*Zt)(Xt_str-C_str*Zt) 
                end
                clearvars t tt
                Q = sumQM/count; % count 1/(t-1), count should correspond to 2 to T
                R_str = sumRM_str/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
                %% %%%%%%%%%%%%%% Test phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                curTrLength = size(stateD{testI},2); % position, velocity in 20ms bins
                % Initialization
                mu_str=Pi;   % sample mean
                sigma_str=V; % sample covariance
                
                clearvars estStateMean_*;
                clearvars estStateCov_*;
                valCellIStr = ~(sum(C_str,2)==0)&~isnan(sum(C_str,2));
                
                C_strVal = C_str(valCellIStr,:); 
                R_strVal = R_str(valCellIStr,valCellIStr); 
                
                for b = 1:curTrLength % # of timebins
                    % One-step prediction by ctx and str data separately
                    mu_str = A*mu_str;
                    sigma_str = A*sigma_str*A'+Q;
                    
                    % compute the Kalman gain (needs to separately computed for ctx and str)
                    K_str = sigma_str*C_strVal'/(C_strVal*sigma_str*C_strVal'+R_strVal); % *inv(C_str*sigma*C_str'+R_str); % Kalman gain for str observation model                 
                    
                    % update the state by str observation
                    curObsX_str = strD{testI}(strI,b); % take the current str spikes bin-by-bin
                    curObsX_str = curObsX_str(valCellIStr,1); 
                    
                    mu_str = mu_str + K_str*(curObsX_str-C_strVal*mu_str);
                    sigma_str = sigma_str - K_str*C_strVal*sigma_str;
                    estStateMean_str(:,b)=mu_str;
                    estStateCov_str(:,:,b)=sigma_str;
                             
                end
                clearvars b
                
                estStateStrMean{r,c,i} = estStateMean_str;
                estStateStrCov{r,c,i} = estStateCov_str;
                             
                %clearvars estStateMean_* estStateCov_*              
            end
        end
    end
    fprintf('finished iteration# %d\n', i)
end
clearvars r c

[strRez.est,strRez.estCut,strRez.estInt] = getRepresentativeStateEst(stateD,estStateStrMean,valTrI);

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
