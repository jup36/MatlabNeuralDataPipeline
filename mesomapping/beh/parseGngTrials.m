function tbytDat = parseGngTrials(tbytDat)

chunkCutoff = 0.4; % 0.4 s licks that occur within this cutoff from one another are chunked together
for jj = 1:length(tbytDat)
    if ~isempty(tbytDat(jj).Lick)
        refTime = tbytDat(jj).stimOn; 

        % Chunk licks
        tbytDat(jj).LickChunk = chunkTimestamps(tbytDat(jj).Lick, chunkCutoff);

        % ITI (find lick chunks where every element belongs to ITI)
        itiChunkI = cell2mat(cellfun(@(a) sum(a < tbytDat(jj).stimOn)==length(a), tbytDat(jj).LickChunk, 'UniformOutput', false));
        tbytDat(jj).itiLickChunk = cellfun(@(a) a-refTime, tbytDat(jj).LickChunk(itiChunkI), 'UniformOutput', false);

        %% classify types of licks in go/no-go trials
        % Hit
        if tbytDat(jj).rewardTrI == 1 && ~isempty(tbytDat(jj).water)
            % find the last lick before water delivery that likely triggered the reward
            lickBeforeWater = tbytDat(jj).Lick(find(tbytDat(jj).Lick < tbytDat(jj).water(1), 1, 'last'));
            hitChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeWater)==1, tbytDat(jj).LickChunk, 'UniformOutput', false));
            hitChunkLicks = tbytDat(jj).LickChunk{hitChunkI};
            tbytDat(jj).hitLicks = hitChunkLicks(hitChunkLicks<=tbytDat(jj).water(1)) - refTime;

            % find licks after the stim offset 
            lickAfterStimOffI = cell2mat(cellfun(@(a) sum(a>=tbytDat(jj).stimOff)==length(a), tbytDat(jj).LickChunk, 'UniformOutput', false)); 
            tbytDat(jj).postStimChunk = cellfun(@(a) a-refTime, tbytDat(jj).LickChunk(lickAfterStimOffI), 'UniformOutput', false);

            % Consumptive licks 
            if sum(tbytDat(jj).Lick > tbytDat(jj).water(1))>0
                % find the first lick after water delivery
                lickAfterWater = tbytDat(jj).Lick(find(tbytDat(jj).Lick > tbytDat(jj).water(1), 1, 'first'));
                consumeChunkI = cell2mat(cellfun(@(a) sum(a==lickAfterWater)==1, tbytDat(jj).LickChunk, 'UniformOutput', false));
                consumeChunkLicks = tbytDat(jj).LickChunk{consumeChunkI};
                tbytDat(jj).consumeLicks = consumeChunkLicks(consumeChunkLicks>tbytDat(jj).water(1)) - refTime;
            end
        % False Alarm
        elseif tbytDat(jj).punishTrI == 1 && ~isempty(tbytDat(jj).airpuff)
            % find the last lick before airpuff that likely triggered the reward
            lickBeforeAir = tbytDat(jj).Lick(find(tbytDat(jj).Lick < tbytDat(jj).airpuff(1), 1, 'last'));
            faChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeAir)==1, tbytDat(jj).LickChunk, 'UniformOutput', false));
            faChunkLicks = tbytDat(jj).LickChunk{faChunkI};
            tbytDat(jj).faLicks = faChunkLicks(faChunkLicks<=tbytDat(jj).airpuff(1)) - refTime;
            
            if sum(tbytDat(jj).Lick > tbytDat(jj).airpuff(1))>0
                % find the first lick after airpuff delivery
                lickAfterAir = tbytDat(jj).Lick(find(tbytDat(jj).Lick > tbytDat(jj).airpuff(1), 1, 'first'));
                postAirChunkI = cell2mat(cellfun(@(a) sum(a==lickAfterAir)==1, tbytDat(jj).LickChunk, 'UniformOutput', false));
                postAirChunkLicks = tbytDat(jj).LickChunk{postAirChunkI};
                tbytDat(jj).postAirLicks = postAirChunkLicks(postAirChunkLicks>tbytDat(jj).airpuff(1)) - refTime;
            end
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chunks = chunkTimestamps(timestamps, cutoff)
% chunk timestamps whose intervals are shorter than the cutoff
    timestamps = timestamps(:); 
    
    % Calculate the difference between consecutive timestamps
    intervals = [Inf; diff(timestamps)];

    % Initialize variables
    chunks = {};
    startIdx = 1; % Start index of the current chunk

    % Iterate through the intervals
    for i = 2:length(timestamps)
        if intervals(i) > cutoff
            % If interval exceeds cutoff, end current chunk and start a new one
            chunks{end + 1} = timestamps(startIdx:i-1);
            startIdx = i;
        end
    end
    
    % Add the last chunk
    chunks{end + 1} = timestamps(startIdx:end);
end



end