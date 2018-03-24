function [ crtseries ] = crttimeseries( rawts, startflg, numbdata )
%This function is to correct the time series data, so that the length of
% the time series (crtts) match the real time (in msec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input                                                                   %   
%   rawts: time series before correction                                  %  
%   startflg: all start timestamps                                        %  
%   numbdata: the number of data points in each chunck of the time series %
% Output                                                                  %
%   crtseries: corrected time series                                      %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

startflag = 1000*(startflg - startflg(1))+1;        % in msec (*1000)
stopflag = startflag + numbdata - 1;

crtseries = nan(1,round(max(stopflag)));

for i = 1:length(startflag)
    if i == 1   % the first fragment 
        crtseries(1, 1:numbdata(1)) = rawts(1, 1:numbdata(1)); 
    elseif i < length(startflg)     % the middle fragment
        crtseries(1, startflag(i):startflag(i)+numbdata(i)-1) = rawts(1,sum(numbdata(1:i-1))+1:sum(numbdata(1:i))); 
    else    % the last fragment
        crtseries(1, startflag(i):startflag(i)+numbdata(i)-1) = rawts(1,sum(numbdata(1:i-1))+1:sum(numbdata(1:i))); 
    end
end

