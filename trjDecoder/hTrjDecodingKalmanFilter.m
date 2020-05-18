
filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles';
cd(filePath)

fileName = 'preprocessKFdecodeHTrjCtxStr_WR40_081919.mat';
load(fullfile(filePath,fileName),'s')

resample = 100;
valTrI = cell2mat(cellfun(@(a) ~isempty(a), s.dat.spkCtx, 'un', 0));
trNumb = sum(valTrI); % # of trials per position/torque pairs
trNumbTot = sum(trNumb); % total number of trial
trainTrN = min(sum(valTrI))-1;
ctxCnumb = max(unique(cell2mat(cellfun(@(a) size(a,1), s.dat.spkCtx, 'un', 0)))); % # of cortex cells
strCnumb = max(unique(cell2mat(cellfun(@(a) size(a,1), s.dat.spkStr, 'un', 0)))); % # of striatal cells
minCnumb = min(ctxCnumb, strCnumb); % # of cells to be included

% leave-a-trial-out decoding using Kalman
for i = 1:resample % repeat resampling trials
    randCtxI = randperm(ctxCnumb);
    randStrI = randperm(strCnumb);
    ctxI = randCtxI(1:minCnumb); % ctx cells for this iteration
    strI = randStrI(1:minCnumb); % str cells for this iteration
    
    for r = 1:size(s.dat.spkCtx,1) % row: trials
        for c = 1:size(s.dat.spkCtx,2) % column: position/torque pairs
            trainI = true(size(valTrI,1),size(valTrI,2));
            if valTrI(r,c)
                trainI(r,c) = false; % to leave one trial out as a test trial
                testI = ~trainI; % index for the one test trial left out
                valTrainI = trainI & valTrI; % valid train trial index
                
                % get current train trials by resampling (to include the same # of trials for each trial type)
                currTrainTrs = zeros(trainTrN,size(trainI,2));
                currTrainState = []; % current training data states (hand position: row 1-3, velocity: row 4-6)
                currTrainCtx = [];
                currTrainStr = [];
                
                for cc = 1:size(valTrainI,2) % sample trials from each trial type
                    tempValTr = find(valTrainI(:,cc));
                    tempValTrRand = tempValTr(randperm(length(tempValTr)));
                    currTrainTrs(:,cc) = tempValTrRand(1:trainTrN); % take the set # of randomized trials from each trial type (column)
                    currTrainState = [currTrainState, s.dat.state(currTrainTrs(:,cc),cc)]; % concatanate randomly selected trials from each type to construct the train data matrix
                    currTrainCtx = [currTrainCtx, cellfun(@(a) a(ctxI,:), s.dat.spkCtx(currTrainTrs(:,cc),cc), 'un', 0)]; % ctx spike mat with cells resampled to match # of cells
                    currTrainStr = [currTrainStr, cellfun(@(a) a(strI,:), s.dat.spkStr(currTrainTrs(:,cc),cc), 'un', 0)]; % str spike mat with cells resampled to match # of cells
                end
                clearvars cc
                [NTrain,NTrType] = size(currTrainState);
                
                %% %%%%%%%%%%% Training phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                interStateM = 0;    % inter-state matrix
                intraStateM = 0;    % intra-state matrix
                obsStateCtxM = 0;   % observed Ctx state matrix
                obsStateStrM = 0;   % observed Str state matrix
                obsIntraStateM = 0; % observed intra-state matrix
                count = 0;
                countPool = 0;
                for t = 1:size(currTrainState,1) % # of trial
                    for tt = 1:size(currTrainState,2) % # of trial types
                        % for parameter A, state model describes how state evolves over time
                        stateZ1 = currTrainState{t,tt}(:,2:end); % Zt
                        stateZ2 = currTrainState{t,tt}(:,1:end-1); % Zt-1
                        tmp1 = stateZ1*stateZ2'; % Zt*Zt-1' (6-by-6 matrix)
                        interStateM = interStateM+tmp1; % sum Zt*Zt-1'
                        tmp2 = stateZ2*stateZ2'; % Zt-1*Zt-1' (6-by-6 matrix)
                        intraStateM = intraStateM+tmp2; % sum Zt-1*Zt-1'
                        % for parameter C, observation model describes how observation relates to the state
                        stateZ = currTrainState{t,tt}; % state (position, velocity)
                        obsCtx = currTrainCtx{t,tt}; % ctx spike data
                        obsStr = currTrainStr{t,tt}; % str spike data
                        tmp3 = obsCtx*stateZ'; % Xt*Zt' (Ctx cell-by-state mat)
                        obsStateCtxM = obsStateCtxM+tmp3; % sum Ctx observation-state matrix
                        tmp4 = obsStr*stateZ'; % Xt*Zt' (Str cell-by-state mat)
                        obsStateStrM = obsStateStrM+tmp4; % sum Str observation-state matrix
                        tmp5 = stateZ*stateZ'; % Zt*Zt'
                        obsIntraStateM = obsIntraStateM+tmp5; % sum Zt*Zt'
                        % for parameter Pi and V
                        count=count+1;
                        poolZStart(:,count) = currTrainState{t,tt}(:,1);
                    end
                end
                clearvars t tt
                A=(interStateM/count)/(intraStateM/count); % Slope for the state model, which defines how state evolves over time
                C_ctx=(obsStateCtxM/count)/(obsIntraStateM/count); % Slope for the observation model for ctx spikes, which defines how observation relates to the state
                C_str=(obsStateStrM/count)/(obsIntraStateM/count); % Slope for the observation model for str spikes
                Pi=mean(poolZStart,2); % Mean for the initial state (sample mean)
                V=cov(poolZStart');    % Covariance for the initial state (sample covariance)
                clearvars poolZStart
                
                % Fit Q and R (variance)
                sumQM = 0;
                sumRM_ctx = 0;
                sumRM_str = 0;
                count = 0;
                for t = 1:size(currTrainState,1) % # of trial
                    for tt = 1:size(currTrainState,2) % # of trial types
                        % for parameter Q, covariance of the state model
                        stateZ1 = currTrainState{t,tt}(:,2:end); % Zt
                        stateZ2 = currTrainState{t,tt}(:,1:end-1); % Zt-1
                        tmp6 = stateZ1-A*stateZ2; % Zt-A*Zt-1
                        tmp6 = tmp6*tmp6'; % (Zt-A*Zt-1)(Zt-A*Zt-1)'
                        count = count + size(stateZ1,2);
                        sumQM = sumQM + tmp6; % Sum (Zt-A*Zt-1)(Zt-A*Zt-1)'
                        % for parameter R of Ctx spikes, covariance of the observation model
                        stateZ = currTrainState{t,tt}; % Zt
                        obsCtx = currTrainCtx{t,tt}; % ctx spike data
                        tmp7 = obsCtx-C_ctx*stateZ; % Xt_ctx-C_ctx*Zt
                        tmp7 = tmp7*tmp7'; % (Xt_ctx-C_ctx*Zt)(Xt_ctx-C_ctx*Zt)
                        sumRM_ctx = sumRM_ctx+tmp7; % Sum (Xt_ctx-C_ctx*Zt)(Xt_ctx-C_ctx*Zt)
                        % for parameter R of Str spikes, covariance of the observation model
                        obsStr = currTrainStr{t,tt}; % str spike data
                        tmp8 = obsStr-C_str*stateZ; % Xt_str-C_str*Zt
                        tmp8 = tmp8*tmp8'; % (Xt_str-C_str*Zt)(Xt_str-C_str*Zt)
                        sumRM_str = sumRM_str+tmp8; % Sum (Xt_str-C_str*Zt)(Xt_str-C_str*Zt)
                    end
                end
                clearvars t tt
                Q = sumQM/count; % count 1/(t-1), count should correspond to 2 to T
                R_ctx = sumRM_ctx/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
                R_str = sumRM_str/(count+NTrain*NTrType); % count 1/t, count should correspond to 1 to T
                
                %% %%%%%%%%%%%%%% Test phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                curTrLength = size(s.dat.state{testI},2); % position, velocity in 20ms bins
                % Initialization
                mu_ctx=Pi;   % sample mean
                sigma_ctx=V; % sample covariance
                mu_str=Pi;   % sample mean
                sigma_str=V; % sample covariance
                clearvars curEstStateMean*;
                clearvars curEstStateCov*;
                for b = 1:curTrLength % # of timebins
                    % One-step prediction by ctx and str data separately
                    mu_ctx = A*mu_ctx; % A is the coefficient matrix for the state model that maps the previous states to current states
                    sigma_ctx = A*sigma_ctx*A'+Q; % Q is the covariance of the state model
                    mu_str = A*mu_str;
                    sigma_str = A*sigma_str*A'+Q;
                    % compute the Kalman gain (needs to separately computed for ctx and str)
                    K_ctx = sigma_ctx*C_ctx'/(C_ctx*sigma_ctx*C_ctx'+R_ctx); % *inv(C_ctx*sigma*C_ctx'+R_ctx); % Kalman gain for ctx observation model
                    K_str = sigma_str*C_str'/(C_str*sigma_str*C_str'+R_str); % *inv(C_str*sigma*C_str'+R_str); % Kalman gain for str observation model
                    % update the state by ctx observation
                    curObsX_ctx = s.dat.spkCtx{testI}(ctxI,b); % take the current ctx spikes bin-by-bin
                    mu_ctx = mu_ctx + K_ctx*(curObsX_ctx-C_ctx*mu_ctx);
                    sigma_ctx = sigma_ctx - K_ctx*C_ctx*sigma_ctx;
                    curEstStateMean_ctx(:,b)=mu_ctx;
                    curEstStateCov_ctx(:,:,b)=sigma_ctx;
                    
                    % update the state by str observation
                    curObsX_str = s.dat.spkStr{testI}(strI,b); % take the current str spikes bin-by-bin
                    mu_str = mu_str + K_str*(curObsX_str-C_str*mu_str);
                    sigma_str = sigma_str - K_str*C_str*sigma_str;
                    curEstStateMean_str(:,b)=mu_str;
                    curEstStateCov_str(:,:,b)=sigma_str;
                end
                s.dat.estStateCtxMean{r,c,i} = curEstStateMean_ctx;
                s.dat.estStateCtxCov{r,c,i} = curEstStateCov_ctx;
                
                s.dat.estStateStrMean{r,c,i} = curEstStateMean_str;
                s.dat.estStateStrCov{r,c,i} = curEstStateCov_str;
                
            end
        end
    end
    fprintf('finished iteration# %d\n', i)
end
clearvars r c i

%save(fullfile(filePath,'rezKFdecodeHTrjCtxStr_WR40_081919.mat'),'s')
%load(fullfile(filePath,'rezKFdecodeHTrjCtxStr_WR40_081919.mat'),'s')
%cell2mat(s.dat.estStateCtxMean(1,1,:))

%% take average across estimated trajectories with resampling
for rr = 1:size(s.dat.estStateCtxMean,1)
    for cc = 1:size(s.dat.estStateCtxMean,2)
        s.dat.avgEstStateCtxMean{rr,cc} = nanmean(cell2mat(s.dat.estStateCtxMean(rr,cc,:)),3); % average across ctx resampled trials
        if size(s.dat.avgEstStateCtxMean{rr,cc},2)>=5
            if size(s.dat.avgEstStateCtxMean{rr,cc},2)>=100
                intAvgEstStateCtxMean{rr,cc} = s.dat.avgEstStateCtxMean{rr,cc}(:,1:100);
            elseif size(s.dat.avgEstStateCtxMean{rr,cc},2)<100
                intAvgEstStateCtxMean{rr,cc} = intm(s.dat.avgEstStateCtxMean{rr,cc},100);
            end
        elseif size(s.dat.avgEstStateCtxMean{rr,cc},2)<5
            intAvgEstStateCtxMean{rr,cc} = nan(size(s.dat.avgEstStateCtxMean{rr,cc},1),100);
        end
        
        s.dat.avgEstStateStrMean{rr,cc} = nanmean(cell2mat(s.dat.estStateStrMean(rr,cc,:)),3); % average across str resampled trials
        if size(s.dat.avgEstStateStrMean{rr,cc},2)>=5
            if size(s.dat.avgEstStateStrMean{rr,cc},2)>=100
                intAvgEstStateStrMean{rr,cc} = s.dat.avgEstStateStrMean{rr,cc}(:,1:100);
            elseif size(s.dat.avgEstStateStrMean{rr,cc},2)<100
                intAvgEstStateStrMean{rr,cc} = intm(s.dat.avgEstStateStrMean{rr,cc},100);
            end
        elseif size(s.dat.avgEstStateStrMean{rr,cc},2)<5
            intAvgEstStateStrMean{rr,cc} = nan(size(s.dat.avgEstStateStrMean{rr,cc},1),100);
        end
        
        if size(s.dat.state{rr,cc},2)>=5
            if size(s.dat.state{rr,cc},2)>=100
                intState{rr,cc} = s.dat.state{rr,cc}(:,1:100);
            elseif size(s.dat.state{rr,cc},2)<100
                intState{rr,cc} = intm(s.dat.state{rr,cc},100);
            end
        elseif size(s.dat.state{rr,cc},2)<5
            intState{rr,cc} = nan(size(s.dat.state{rr,cc},1),100);
        end
        
        % distance from the actual trajectory
        s.dat.dTrjCtx{rr,cc} = sqrt((s.dat.state{rr,cc}-s.dat.avgEstStateCtxMean{rr,cc}).^2);
        s.dat.dTrjStr{rr,cc} = sqrt((s.dat.state{rr,cc}-s.dat.avgEstStateStrMean{rr,cc}).^2);
    end
end
clearvars rr cc

%% Plot interpolated trajectories
% low Torque Position 1 (column 1)
ctxLtP1 = intAvgEstStateCtxMean(:,1);
ctxLtP1 = ctxLtP1(~cellfun(@isempty, ctxLtP1)); 
ctxLtP1X = cell2mat(cellfun(@(a) a(1,:), ctxLtP1(:,1), 'un', 0)); 

strLtP1 = intAvgEstStateStrMean(:,1);
strLtP1 = strLtP1(~cellfun(@isempty, strLtP1)); 
strLtP1X = cell2mat(cellfun(@(a) a(1,:), strLtP1(:,1), 'un', 0)); 

actLtP1 = intState(:,1);
actLtP1 = actLtP1(~cellfun(@isempty, actLtP1)); 
actLtP1X = cell2mat(cellfun(@(a) a(1,:), actLtP1(:,1), 'un', 0)); 

hold on; 
plot(mean(ctxLtP1X),'b'); 
plot(mean(strLtP1X),'r'); 
plot(mean(actLtP1X),'g'); 
hold off; 

% low Torque Position 2 (column 3)
ctxLtP2 = intAvgEstStateCtxMean(:,3);
ctxLtP2 = ctxLtP2(~cellfun(@isempty, ctxLtP2)); 
ctxLtP2X = cell2mat(cellfun(@(a) a(1,:), ctxLtP2(:,1), 'un', 0)); 

strLtP2 = intAvgEstStateStrMean(:,3);
strLtP2 = strLtP2(~cellfun(@isempty, strLtP2)); 
strLtP2X = cell2mat(cellfun(@(a) a(1,:), strLtP2(:,1), 'un', 0)); 

actLtP2 = intState(:,3);
actLtP2 = actLtP2(~cellfun(@isempty, actLtP2)); 
actLtP2X = cell2mat(cellfun(@(a) a(1,:), actLtP2(:,1), 'un', 0)); 

hold on; 
plot(mean(ctxLtP2X),'b'); 
plot(mean(strLtP2X),'r'); 
plot(mean(actLtP2X),'g'); 
hold off; 

% high Torque Position 1 (column 2)
ctxHtP1 = intAvgEstStateCtxMean(:,2);
ctxHtP1 = ctxHtP1(~cellfun(@isempty, ctxHtP1)); 
ctxHtP1X = cell2mat(cellfun(@(a) a(1,:), ctxHtP1(:,1), 'un', 0)); 

strHtP1 = intAvgEstStateStrMean(:,2);
strHtP1 = strHtP1(~cellfun(@isempty, strHtP1)); 
strHtP1X = cell2mat(cellfun(@(a) a(1,:), strHtP1(:,1), 'un', 0)); 

actHtP1 = intState(:,2);
actHtP1 = actHtP1(~cellfun(@isempty, actHtP1)); 
actHtP1X = cell2mat(cellfun(@(a) a(1,:), actHtP1(:,1), 'un', 0)); 

hold on; 
plot(mean(ctxHtP1X),'b'); 
plot(mean(strHtP1X),'r'); 
plot(mean(actHtP1X),'g'); 
hold off; 

% high Torque Position 2 (column 4)
ctxHtP2 = intAvgEstStateCtxMean(:,4);
ctxHtP2 = ctxHtP2(~cellfun(@isempty, ctxHtP2)); 
ctxHtP2X = cell2mat(cellfun(@(a) a(1,:), ctxHtP2(:,1), 'un', 0)); 

strHtP2 = intAvgEstStateStrMean(:,4);
strHtP2 = strHtP2(~cellfun(@isempty, strHtP2)); 
strHtP2X = cell2mat(cellfun(@(a) a(1,:), strHtP2(:,1), 'un', 0)); 

actHtP2 = intState(:,4);
actHtP2 = actHtP2(~cellfun(@isempty, actHtP2)); 
actHtP2X = cell2mat(cellfun(@(a) a(1,:), actHtP2(:,1), 'un', 0)); 

hold on; 
plot(mean(ctxHtP2X),'b'); 
plot(mean(strHtP2X),'r'); 
plot(mean(actHtP2X),'g'); 
hold off; 

% meanLtP1vsP2_CtxStrActual
hold on;
plot(mean(actLtP1X),'g');
plot(mean(actLtP2X),'g');
plot(mean(ctxLtP1X),'b');
plot(mean(ctxLtP2X),'b');
plot(mean(strLtP1X),'r');
plot(mean(strLtP2X),'r');
hold off; 
print(fullfile(filePath,'Figure','meanLtP1vsP2_CtxStrActual.pdf'),'-bestfit','-dpdf','-painters')

% dmeanLtP1vsP2_CtxStrActual
hold on;
plot(mean(actLtP1X)-mean(actLtP2X),'g');
plot(mean(ctxLtP1X)-mean(ctxLtP2X),'b');
plot(mean(strLtP1X)-mean(strLtP2X),'r');
hold off; 
xlim([1 50])
print(fullfile(filePath,'Figure','dmeanLtP1vsP2_CtxStrActual.pdf'),'-bestfit','-dpdf','-painters')


%% plot concatenated
timeX = 0;
figure; hold on;
for c = 1:size(s.dat.estStateCtxMean,2)
    for r = 1:20 % just to include the first block only per trial type
        if ~isempty(s.dat.state{r,c})
            ctxTrjX = s.dat.avgEstStateCtxMean{r,c}(1,:); % Ctx X trj (left-right, horizontal hand position)
            strTrjX = s.dat.avgEstStateStrMean{r,c}(1,:); % Str X trj
            actTrjX = s.dat.state{r,c}(1,:); % actual trj
            
            tempX = timeX+2:timeX+length(actTrjX)+1;
            timeX = timeX+length(actTrjX)+1; % update timeX
            
            if s.dat.pos1{r,c}==min([s.dat.pos1{:}]) % left
                if s.dat.trq{r,c}==min([s.dat.trq{:}]) % low-torque
                    plot(tempX,actTrjX,'k','LineWidth',1);
                    plot(tempX,ctxTrjX,'Color',[100 149 237]./255,'LineWidth',1);
                    plot(tempX,strTrjX,'Color',[50 205 50]./255,'LineWidth',1);
                elseif s.dat.trq{r,c}==max([s.dat.trq{:}]) % high-torque
                    plot(tempX,actTrjX,'k','LineWidth',2);
                    plot(tempX,ctxTrjX,'Color',[100 149 237]./255,'LineWidth',2);
                    plot(tempX,strTrjX,'Color',[50 205 50]./255,'LineWidth',2);
                end
            elseif s.dat.pos1{r,c}==max([s.dat.pos1{:}]) % right
                if s.dat.trq{r,c}==min([s.dat.trq{:}]) % low-torque
                    plot(tempX,actTrjX,'k','LineWidth',1);
                    plot(tempX,ctxTrjX,'Color',[100 149 237]./255,'LineWidth',1);
                    plot(tempX,strTrjX,'Color',[50 205 50]./255,'LineWidth',1);
                elseif s.dat.trq{r,c}==max([s.dat.trq{:}]) % high-torque
                    plot(tempX,actTrjX,'k','LineWidth',2);
                    plot(tempX,ctxTrjX,'Color',[100 149 237]./255,'LineWidth',2);
                    plot(tempX,strTrjX,'Color',[50 205 50]./255,'LineWidth',2);
                end
            end
        end
    end
end
hold off;
clearvars r c

%% compute residuals errors (1. from the beginning, 2. after pull stop)
for r = 1:size(s.dat.estStateCtxMean,1)
    for c = 1:size(s.dat.estStateCtxMean,2)
        if ~isempty(s.dat.state{r,c})
            % get pullStop point
            pStop = find(s.dat.pullIdx{r,c},1,'last');
            if isempty(pStop)
                pStop = size(s.dat.pullIdx{r,c},2);
            end
            rsdCtx1{r,c} = s.dat.avgEstStateCtxMean{r,c}(:,1:pStop)-s.dat.state{r,c}(:,1:pStop);
            rsdCtx1{r,c} = sqrt(rsdCtx1{r,c}.^2); % squared error
            % interpolation
            if size(rsdCtx1{r,c},2)>=5
                if size(rsdCtx1{r,c},2)>=100
                    rsdCtx1Int{r,c} = rsdCtx1{r,c}(:,1:100);
                elseif size(rsdCtx1{r,c},2)<100
                    rsdCtx1Int{r,c} = intm(rsdCtx1{r,c},100);
                end
            elseif size(rsdCtx1{r,c},2)<5
                rsdCtx1Int{r,c} = nan(size(rsdCtx1{r,c},1),100);
            end
            
            % interpolation
            rsdCtx2{r,c} = s.dat.avgEstStateCtxMean{r,c}(:,pStop:end)-s.dat.state{r,c}(:,pStop:end);
            rsdCtx2{r,c} = sqrt(rsdCtx2{r,c}.^2);
            if size(rsdCtx2{r,c},2)>=5
                if size(rsdCtx2{r,c},2)>=100
                    rsdCtx2Int{r,c} = rsdCtx2{r,c}(:,1:100);
                elseif size(rsdCtx2{r,c},2)<100
                    rsdCtx2Int{r,c} = intm(rsdCtx2{r,c},100);
                end
            elseif size(rsdCtx2{r,c},2)<5
                rsdCtx2Int{r,c} = nan(size(rsdCtx2{r,c},1),100);
            end
            
            % interpolation
            rsdStr1{r,c} = s.dat.avgEstStateStrMean{r,c}(:,1:pStop)-s.dat.state{r,c}(:,1:pStop);
            rsdStr1{r,c} = sqrt(rsdStr1{r,c}.^2);
            if size(rsdStr1{r,c},2)>=5
                if size(rsdStr1{r,c},2)>=100
                    rsdStr1Int{r,c} = rsdStr1{r,c}(:,1:100);
                elseif size(rsdStr1{r,c},2)<100
                    rsdStr1Int{r,c} = intm(rsdStr1{r,c},100);
                end
            elseif size(rsdStr1{r,c},2)<5
                rsdStr1Int{r,c} = nan(size(rsdStr1{r,c},1),100);
            end
            
            % interpolation
            rsdStr2{r,c} = s.dat.avgEstStateStrMean{r,c}(:,pStop:end)-s.dat.state{r,c}(:,pStop:end);
            rsdStr2{r,c} = sqrt(rsdStr2{r,c}.^2);
            
            if size(rsdStr2{r,c},2)>=5
                if size(rsdStr2{r,c},2)>=100
                    rsdStr2Int{r,c} = rsdStr2{r,c}(:,1:100);
                elseif size(rsdStr2{r,c},2)<100
                    rsdStr2Int{r,c} = intm(rsdStr2{r,c},100);
                end
            elseif size(rsdStr2{r,c},2)<5
                rsdStr2Int{r,c} = nan(size(rsdStr2{r,c},1),100);
            end
        end
    end
end
clearvars r c

msRsdToPullStopCtx = cell(2,size(s.dat.state,2));  
msRsdFrPullStopCtx = cell(2,size(s.dat.state,2));  
msRsdToPullStopStr = cell(2,size(s.dat.state,2));  
msRsdFrPullStopStr = cell(2,size(s.dat.state,2));  

for c = 1:size(rsdCtx1Int,2)
    
    currToPullStopCtx = reshape([rsdCtx1Int{:,c}],size(rsdCtx1Int{c},1),size(rsdCtx1Int{c},2),[]); 
    currFrPullStopCtx = reshape([rsdCtx2Int{:,c}],size(rsdCtx2Int{c},1),size(rsdCtx2Int{c},2),[]); 
    
    currToPullStopStr = reshape([rsdStr1Int{:,1}],size(rsdStr1Int{1},1),size(rsdStr1Int{1},2),[]); 
    currFrPullStopStr = reshape([rsdStr2Int{:,1}],size(rsdStr2Int{1},1),size(rsdStr2Int{1},2),[]); 
    
    for tt = 1:6 % increament kinematic variables (x,y,z) x 2 (position, velocity)
        [msRsdToPullStopCtx{1,c}(tt,:),~,msRsdToPullStopCtx{2,c}(tt,:)] = meanstdsem(squeeze(currToPullStopCtx(tt,:,:))');
        [msRsdFrPullStopCtx{1,c}(tt,:),~,msRsdFrPullStopCtx{2,c}(tt,:)] = meanstdsem(squeeze(currFrPullStopCtx(tt,:,:))');  
        
        [msRsdToPullStopStr{1,c}(tt,:),~,msRsdToPullStopStr{2,c}(tt,:)] = meanstdsem(squeeze(currToPullStopStr(tt,:,:))');
        [msRsdFrPullStopStr{1,c}(tt,:),~,msRsdFrPullStopStr{2,c}(tt,:)] = meanstdsem(squeeze(currFrPullStopStr(tt,:,:))');  
    end
end 
clearvars c tt

%% plot residual
colorMap = [[100 149 237]./255; [50 205 50]./255];
figure; boundedline(1:100, msRsdToPullStopCtx{1,1}(1,:), msRsdToPullStopCtx{2,1}(1,:),...
    1:100, msRsdToPullStopStr{1,1}(1,:), msRsdToPullStopStr{2,1}(1,:), 'cmap', colorMap, 'transparency', 0.2);

figure; boundedline(1:100, msRsdToPullStopCtx{1,1}(2,:), msRsdToPullStopCtx{2,1}(2,:),...
    1:100, msRsdToPullStopStr{1,1}(2,:), msRsdToPullStopStr{2,1}(1,:), 'cmap', colorMap, 'transparency', 0.2);

figure; boundedline(1:100, msRsdToPullStopCtx{1,1}(3,:), msRsdToPullStopCtx{2,1}(3,:),...
    1:100, msRsdToPullStopStr{1,1}(3,:), msRsdToPullStopStr{2,1}(3,:), 'cmap', colorMap, 'transparency', 0.2);

figure; boundedline(1:100, msRsdToPullStopCtx{1,1}(4,:), msRsdToPullStopCtx{2,1}(4,:),...
    1:100, msRsdToPullStopStr{1,1}(4,:), msRsdToPullStopStr{2,1}(4,:), 'cmap', colorMap, 'transparency', 0.2);

figure; boundedline(1:100, msRsdToPullStopCtx{1,1}(5,:), msRsdToPullStopCtx{2,1}(5,:),...
    1:100, msRsdToPullStopStr{1,1}(5,:), msRsdToPullStopStr{2,1}(5,:), 'cmap', colorMap, 'transparency', 0.2);

figure; boundedline(1:100, msRsdToPullStopCtx{1,1}(6,:), msRsdToPullStopCtx{2,1}(6,:),...
    1:100, msRsdToPullStopStr{1,1}(6,:), msRsdToPullStopStr{2,1}(6,:), 'cmap', colorMap, 'transparency', 0.2);





% Matrix interpolation function 
function [intMat] = intm(origMat, numbDataPoints)
% numbDataPoints = 100; % 20ms*100 = 2000ms
if numbDataPoints>size(origMat,2)
    x=1:size(origMat,2);
    xq=linspace(1,size(origMat,2),numbDataPoints);
    intMat = interp1(x,origMat',xq)';
else
    intMat = origMat;
end
end






















