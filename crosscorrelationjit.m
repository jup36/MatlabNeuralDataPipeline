function [ xcind, xcsum, confbnd, lag ] = crosscorrelationjit( spktrn1, spktrn2, uv )
%crosscorrelation performs crosscorrelation analysis using xcorr
% Input: spike train pair 

for t = 1:uv.trnumb % increment trials  
    xc(t,:) = xcorr(spktrn1(t,:), spktrn2(t,:), uv.maxlag); % get the cross-correlogram
end
clearvars t temp* tp*
xcsum = sum(xc,1);

shufxc = zeros(100,size(xc,2));

shuf = [-10:10]; % jitter factor -10 to 10 ms
for s = 1:100 % 100 times
    rdxc = zeros(uv.trnumb,size(xcsum,2));
    for t = 1:uv.trnumb %uv.trnumb      
        % get the jittered spike train for spktrn1
        spks1 = find(spktrn1(t,:)==1); % current binned pfc spike train
        shufspks1 = datasample(shuf,length(spks1)) + spks1; % jitter the spikes
        
        while sum(shufspks1 <= 0)>0 || sum(shufspks1 > length(spktrn1))>0 % repeat if jittering got a spike out of the range
              shufspks1 = datasample(shuf,length(spks1)) + spks1;
        end
        
        rdspktrn1 = zeros(1,length(spktrn1));
        rdspktrn1(shufspks1) = 1;
        
        % get the jittered spike train for spktrn2
        spks2 = find(spktrn2(t,:)==1); % current binned pfc spike train
        shufspks2 = datasample(shuf,length(spks2)) + spks2;
        
        while sum(shufspks2 <= 0)>0 || sum(shufspks2 >= length(spktrn2))>0
              shufspks2 = datasample(shuf,length(spks2)) + spks2;
        end
        
        rdspktrn2 = zeros(1,length(spktrn2));
        rdspktrn2(shufspks2) = 1;     
        
        % 
        rdxc(t,:) = xcorr(rdspktrn1,rdspktrn2,uv.maxlag); % get the cross-correlogram with the jittered spike trains
 end
    shufxc(s,:)=sum(rdxc,1);
end
clearvars t s

% find the 99% confidence interval
confbnd = zeros(4,size(shufxc,2));
confbnd(1,:) = max(mean(shufxc,1) + 2.58*std(shufxc,0,1)); % global-wise 99% interval
confbnd(2,:) = mean(shufxc,1) + 2.58*std(shufxc,0,1); % point-wise upper 99% confidence interval
confbnd(3,:) = mean(shufxc,1) - 2.58*std(shufxc,0,1); % point-wise bottom 99% confidence interval
confbnd(4,:) = min(mean(shufxc,1) - 2.58*std(shufxc,0,1)); % global-wise 99% interval

% find the time lags crossing the cofidence interval
timelag = [-uv.maxlag:1:uv.maxlag];
if sum(xcsum > confbnd(1,:)) > 0 && sum(xcsum < confbnd(4,:)) > 0
    xcind = 2;  % crosses both bounds
    lag = [timelag(xcsum < confbnd(4,:)), timelag(xcsum > confbnd(1,:))]; % crossing timelags
elseif sum(xcsum > confbnd(1,:)) > 0
    xcind = 1;  % crosses upperbound
    lag = timelag(xcsum > confbnd(1,:)); % crossing timelags
elseif sum(xcsum < confbnd(4,:)) > 0
    xcind = -1; % crosses bottombound
    lag = timelag(xcsum < confbnd(4,:)); % crossing timelags
else % no crossing
    xcind = 0;
    lag = NaN;
end
     

end

