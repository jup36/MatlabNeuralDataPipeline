    function zscoreMat = binTimestampsZscoreEachTertile(timeStampC, var)

        % Example lick timestamps for 100 trials
        %numTrials = 100;
        %lick_timestamps_trials = cell(1, numTrials);
        %for i = 1:numTrials
        %   lick_timestamps_trials{i} = rand(1, 100)*15 - 5; % random licks between -5s to 10s
        %end

        numBlocks = length(timeStampC);

        for bl = 1:numBlocks
            zscoreC{bl} = binTimestampsZscore(timeStampC{bl}, var);
        end

        zscoreMat = cell2mat([zscoreC(:)]);

    end
