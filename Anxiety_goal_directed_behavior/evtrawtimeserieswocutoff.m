%This function 'evtrawtimeseries' is to get the raw LFP timeseries time-locked 
% to task events.

function [ current_rawmat, current_rawbasemat ] = evtrawtimeserieswocutoff( currentd, currenttstp, currentcue, baseendpoint )

if isnan(currenttstp) == 0      % in case the current tstp is valid (non-NaN)
    
    tempd = createdatamatc(currentd,currenttstp,1000,[2.25 2.25]);      % get the temporary data using the createdatamatc, +- .25 are to give extra spacing due to the moving window (if the size of moving windows is 500 ms, start at -250 ms becuase the mtspecgrams takes 500 ms from the starting point).
    tempcue = currentcue(max(find(currentcue <= currenttstp)));               % find the ITI of each trial
    tempbased = createdatamatc(currentd,tempcue+baseendpoint,1000,[4.5 0]);  % get the temporary base data using the createdatamatc, sv.baseendpoint is set to 9.25 sec considering the moving window
       
    %% this is the new criterion for baseline, in which the entire baseline will be thrown out if there's any violation of the upper or lower limit.
    current_rawmat = tempd';
    current_rawbasemat = tempbased';

else        % in case the current tstp is NaN
    current_rawmat = nan(1,length(tempd))';
    current_rawbasemat = nan(1,length(tempd))';
end

end



