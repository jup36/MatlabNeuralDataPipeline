function [timeStampOut, expd, intTime] = lickSim(timeStampC, minMaxTime, numReps)
        %This function takes a cell array comprising timestamps. The
        % timestamps are turned into interval times to be fitted with an
        % exponential distribution. The fitted exponential pdf is used to generate
        % interval times that are turned back into the timestamps. 
        % timeStampC: a cell array that contains timestamps of which interval
        % times to be modeled. 
        % timeStampC = lickTimeEarlyCue;
        % minMaxTime = [0 5];
        % minTime = 0; % sec
        % numReps = 1000; 

        % sanity check: ensure that all timestamps are within the min and max range
        timeStampC = cellfun(@(a) a(a >= minMaxTime(1) & a <= minMaxTime(2)), timeStampC, 'UniformOutput', false);

        intTimeC = cellfun(@(a) diff([minMaxTime(1); a]), timeStampC, 'UniformOutput', false); % transform to interval time
        intTime = cell2mat(intTimeC'); % take all interval times 

        % fit exponential distribution
        expd = fitdist(intTime, 'Exponential');

        lambda = 1/mean(intTime); % lambda is the rate parameter of the exponential distribution

        timeStampOut = cell(1, numReps); 
        % generate timestamps
        for jj = 1:numReps
            timeStampOut{jj} = expTimestampsTrial(expd, minMaxTime); 
        end


        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function timeStamps = expTimestampsTrial(expd, minMaxTime)
        minTime = min(minMaxTime); 
        maxTime = max(minMaxTime); 
     
        curTime = minTime; % starting point
        intvs = []; % generated intervals
        while curTime <= maxTime
            intv = random(expd, 1, 1);
            curTime = curTime + intv; 
            if curTime <= maxTime
                intvs = [intvs; intv]; 
            end
        end
        % generate cumulative timestamps
        timeStamps = cumsum(intvs) + minTime; 
        end
        
    end