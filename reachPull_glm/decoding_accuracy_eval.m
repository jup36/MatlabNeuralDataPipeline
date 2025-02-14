isStr = cell2mat(spkTimesCell(5,:)); 

%% decoding accuracy - both cortex, striatum
for tt = 1:length(posTqTypes) % trial types
    % log sum probability (target position) 
    logSum_probLe_ext = squeeze(sum(log(cell2mat(decode.probLe(tt,:,:))),3));   % logSum probability of left
    logSum_probRi_ext = squeeze(sum(log(1-cell2mat(decode.probLe(tt,:,:))),3)); % logSum probability of right
    
    logSum_probLe = logSum_probLe_ext(in_t(401:600),:);
    logSum_probRi = logSum_probRi_ext(in_t(401:600),:);
    
    if contains(posTqTypes{tt},'p1')
        decode.targPosAccuracy{tt} = mean(logSum_probLe>=logSum_probRi,2); % correct if left
    elseif contains(posTqTypes{tt},'p2')
        decode.targPosAccuracy{tt} = mean(logSum_probLe<=logSum_probRi,2); % correct if right
    end
    
    % log sum probability (torque) 
    logSum_probLo_ext = squeeze(sum(log(cell2mat(decode.probLo(tt,:,:))),3));   % logSum probability of low
    logSum_probHi_ext = squeeze(sum(log(1-cell2mat(decode.probLo(tt,:,:))),3)); % logSum probability of high
    
    logSum_probLo = logSum_probLo_ext(in_t(401:600),:);
    logSum_probHi = logSum_probHi_ext(in_t(401:600),:);
    
    if contains(posTqTypes{tt},'t1')
        decode.torqAccuracy{tt} = mean(logSum_probLo>=logSum_probHi,2); % correct if left
    elseif contains(posTqTypes{tt},'t2')
        decode.torqAccuracy{tt} = mean(logSum_probLo<=logSum_probHi,2); % correct if right
    end
end

% decode reward or not
for tt = 1:2
    % log sum probability (reward or no) 
    logSum_probNoRwd_ext = squeeze(sum(log(cell2mat(decode.probNoRwd(tt,:,:))),3));   % logSum probability of no-reward
    logSum_probRwd_ext   = squeeze(sum(log(1-cell2mat(decode.probNoRwd(tt,:,:))),3)); % logSum probability of reward
    
    logSum_probNoRwd = logSum_probNoRwd_ext(in_t(401:600),:);
    logSum_probRwd   = logSum_probRwd_ext(in_t(401:600),:);
    
    if tt==1
        decode.rwdAccuracy{tt} = mean(logSum_probNoRwd>=logSum_probRwd,2); % correct if no-reward
    elseif tt==2
        decode.rwdAccuracy{tt} = mean(logSum_probNoRwd<=logSum_probRwd,2); % correct if reward
    end  
end

%% decoding accuracy - cortex
for tt = 1:length(posTqTypes) % trial types
    % log sum probability
    logSum_probLe_ext = squeeze(sum(log(cell2mat(decode.probLe(tt,:,~isStr))),3));   % logSum probability of left
    logSum_probRi_ext = squeeze(sum(log(1-cell2mat(decode.probLe(tt,:,~isStr))),3)); % logSum probability of right
    
    logSum_probLe = logSum_probLe_ext(in_t(401:600),:);
    logSum_probRi = logSum_probRi_ext(in_t(401:600),:);
    
    if contains(posTqTypes{tt},'p1')
        decode.targPosAccuracy_ctx{tt} = mean(logSum_probLe>=logSum_probRi,2); % correct if left
    elseif contains(posTqTypes{tt},'p2')
        decode.targPosAccuracy_ctx{tt} = mean(logSum_probLe<=logSum_probRi,2); % correct if right
    end
    
    % log sum probability
    logSum_probLo_ext = squeeze(sum(log(cell2mat(decode.probLo(tt,:,~isStr))),3));   % logSum probability of low
    logSum_probHi_ext = squeeze(sum(log(1-cell2mat(decode.probLo(tt,:,~isStr))),3)); % logSum probability of high
    
    logSum_probLo = logSum_probLo_ext(in_t(401:600),:);
    logSum_probHi = logSum_probHi_ext(in_t(401:600),:);
    
    if contains(posTqTypes{tt},'t1')
        decode.torqAccuracy_ctx{tt} = mean(logSum_probLo>=logSum_probHi,2); % correct if left
    elseif contains(posTqTypes{tt},'t2')
        decode.torqAccuracy_ctx{tt} = mean(logSum_probLo<=logSum_probHi,2); % correct if right
    end  
end

for tt = 1:2
    % log sum probability (reward or no) 
    logSum_probNoRwd_ext = squeeze(sum(log(cell2mat(decode.probNoRwd(tt,:,~isStr))),3));   % logSum probability of no-reward
    logSum_probRwd_ext   = squeeze(sum(log(1-cell2mat(decode.probNoRwd(tt,:,~isStr))),3)); % logSum probability of reward
    
    logSum_probNoRwd = logSum_probNoRwd_ext(in_t(401:600),:);
    logSum_probRwd   = logSum_probRwd_ext(in_t(401:600),:);
    
    if tt==1
        decode.rwdAccuracy_ctx{tt} = mean(logSum_probNoRwd>=logSum_probRwd,2); % correct if no-reward
    elseif tt==2
        decode.rwdAccuracy_ctx{tt} = mean(logSum_probNoRwd<=logSum_probRwd,2); % correct if reward
    end  
end

%% decoding accuracy - striatum
for tt = 1:length(posTqTypes) % trial types
    % log sum probability
    logSum_probLe_ext = squeeze(sum(log(cell2mat(decode.probLe(tt,:,isStr))),3));   % logSum probability of left
    logSum_probRi_ext = squeeze(sum(log(1-cell2mat(decode.probLe(tt,:,isStr))),3)); % logSum probability of right
    
    logSum_probLe = logSum_probLe_ext(in_t(401:600),:);
    logSum_probRi = logSum_probRi_ext(in_t(401:600),:);
    
    if contains(posTqTypes{tt},'p1')
        decode.targPosAccuracy_str{tt} = mean(logSum_probLe>=logSum_probRi,2); % correct if left
    elseif contains(posTqTypes{tt},'p2')
        decode.targPosAccuracy_str{tt} = mean(logSum_probLe<=logSum_probRi,2); % correct if right
    end
    
    % log sum probability
    logSum_probLo_ext = squeeze(sum(log(cell2mat(decode.probLo(tt,:,isStr))),3));   % logSum probability of low
    logSum_probHi_ext = squeeze(sum(log(1-cell2mat(decode.probLo(tt,:,isStr))),3)); % logSum probability of high
    
    logSum_probLo = logSum_probLo_ext(in_t(401:600),:);
    logSum_probHi = logSum_probHi_ext(in_t(401:600),:);
    
    if contains(posTqTypes{tt},'t1')
        decode.torqAccuracy_str{tt} = mean(logSum_probLo>=logSum_probHi,2); % correct if left
    elseif contains(posTqTypes{tt},'t2')
        decode.torqAccuracy_str{tt} = mean(logSum_probLo<=logSum_probHi,2); % correct if right
    end  
end

for tt = 1:2
    % log sum probability (reward or no) 
    logSum_probNoRwd_ext = squeeze(sum(log(cell2mat(decode.probNoRwd(tt,:,isStr))),3));   % logSum probability of no-reward
    logSum_probRwd_ext   = squeeze(sum(log(1-cell2mat(decode.probNoRwd(tt,:,isStr))),3)); % logSum probability of reward
    
    logSum_probNoRwd = logSum_probNoRwd_ext(in_t(401:600),:);
    logSum_probRwd   = logSum_probRwd_ext(in_t(401:600),:);
    
    if tt==1
        decode.rwdAccuracy_str{tt} = mean(logSum_probNoRwd>=logSum_probRwd,2); % correct if no-reward
    elseif tt==2
        decode.rwdAccuracy_str{tt} = mean(logSum_probNoRwd<=logSum_probRwd,2); % correct if reward
    end  
end

%% summary
figure; hold on; 
plot(smooth2a(nanmean(cell2mat(decode.torqAccuracy),2),2,0))
plot(smooth2a(nanmean(cell2mat(decode.torqAccuracy_ctx),2),2,0))
plot(smooth2a(nanmean(cell2mat(decode.torqAccuracy_str),2),2,0))

figure; hold on; 
plot(smooth2a(nanmean(cell2mat(decode.targPosAccuracy),2),2,0))
plot(smooth2a(nanmean(cell2mat(decode.targPosAccuracy_ctx),2),2,0))
plot(smooth2a(nanmean(cell2mat(decode.targPosAccuracy_str),2),2,0))

figure; hold on; 
plot(smooth2a(nanmean(cell2mat(decode.rwdAccuracy),2),2,0))
plot(smooth2a(nanmean(cell2mat(decode.rwdAccuracy_ctx),2),2,0))
plot(smooth2a(nanmean(cell2mat(decode.rwdAccuracy_str),2),2,0))


