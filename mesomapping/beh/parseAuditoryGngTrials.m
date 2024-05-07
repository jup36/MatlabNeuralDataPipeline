function tbytDat = parseAuditoryGngTrials(tbytDat)

tbytDat = tbytDatRewardPunishI(tbytDat);

chunkCutoff = 0.4; % 0.4 s licks that occur within this cutoff from one another are chunked together
for jj = 1:length(tbytDat)
    if ~isempty(tbytDat(jj).Lick)
        refTime = tbytDat(jj).evtOn;

        % Chunk licks
        tbytDat(jj).LickChunk = chunkTimestamps(tbytDat(jj).Lick, chunkCutoff);

        % ITI (find lick chunks where every element belongs to ITI)
        itiChunkI = cell2mat(cellfun(@(a) sum(a < tbytDat(jj).evtOn)==length(a), tbytDat(jj).LickChunk, 'UniformOutput', false));
        tbytDat(jj).itiLickChunk = cellfun(@(a) a-refTime, tbytDat(jj).LickChunk(itiChunkI), 'UniformOutput', false);

        %% classify types of licks in go/no-go trials
        % Hit
        if (tbytDat(jj).pos_rwd_tr == 1 && ~isempty(tbytDat(jj).water)) || ...
                (tbytDat(jj).pos_oms_tr == 1 && sum(tbytDat(jj).Lick>tbytDat(jj).evtOff)>0) || ...
                (tbytDat(jj).pos_pns_tr == 1 && sum(tbytDat(jj).Lick>tbytDat(jj).evtOff)>0)
            % find the last lick before water delivery that likely triggered the reward
            if ~isempty(tbytDat(jj).water)
                lickBeforeWater = tbytDat(jj).Lick(find(tbytDat(jj).Lick < tbytDat(jj).water(1), 1, 'last'));
                hitChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeWater)==1, tbytDat(jj).LickChunk, 'UniformOutput', false));
                hitChunkLicks = tbytDat(jj).LickChunk{hitChunkI};
                tbytDat(jj).hitLicks = hitChunkLicks(hitChunkLicks<=tbytDat(jj).water(1)) - refTime;
            else
                lickBeforeWater = tbytDat(jj).Lick(find(tbytDat(jj).Lick > tbytDat(jj).evtOff, 1, 'first'));
                hitChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeWater)==1, tbytDat(jj).LickChunk, 'UniformOutput', false));
                hitChunkLicks = tbytDat(jj).LickChunk{hitChunkI};
                tbytDat(jj).hitLicks = hitChunkLicks - refTime;
            end

            % find licks after the stim offset
            lickAfterevtOffI = cell2mat(cellfun(@(a) sum(a>=tbytDat(jj).evtOff)==length(a), tbytDat(jj).LickChunk, 'UniformOutput', false));
            tbytDat(jj).postStimChunk = cellfun(@(a) a-refTime, tbytDat(jj).LickChunk(lickAfterevtOffI), 'UniformOutput', false);

            if ~isempty(tbytDat(jj).water)
                % Consumptive licks
                if sum(tbytDat(jj).Lick > tbytDat(jj).water(1))>0
                    % find the first lick after water delivery
                    lickAfterWater = tbytDat(jj).Lick(find(tbytDat(jj).Lick > tbytDat(jj).water(1), 1, 'first'));
                    consumeChunkI = cell2mat(cellfun(@(a) sum(a==lickAfterWater)==1, tbytDat(jj).LickChunk, 'UniformOutput', false));
                    consumeChunkLicks = tbytDat(jj).LickChunk{consumeChunkI};
                    tbytDat(jj).consumeLicks = consumeChunkLicks(consumeChunkLicks>tbytDat(jj).water(1)) - refTime;
                end
            else
                tbytDat(jj).consumeLicks = []; 
            end
            % False Alarm
        elseif (tbytDat(jj).neg_rwd_tr == 1 && ~isempty(tbytDat(jj).water)) || ...
                (tbytDat(jj).neg_oms_tr == 1 && sum(tbytDat(jj).Lick>tbytDat(jj).evtOff)>0) || ...
                (tbytDat(jj).neg_pns_tr == 1 && sum(tbytDat(jj).Lick>tbytDat(jj).evtOff)>0)
            % find the last lick before airpuff that likely triggered the
            % airpuff
            if ~isempty(tbytDat(jj).airpuff)
                lickBeforeAir = tbytDat(jj).Lick(find(tbytDat(jj).Lick < tbytDat(jj).airpuff(1), 1, 'last'));
                faChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeAir)==1, tbytDat(jj).LickChunk, 'UniformOutput', false));
                faChunkLicks = tbytDat(jj).LickChunk{faChunkI};
                tbytDat(jj).faLicks = faChunkLicks(faChunkLicks<=tbytDat(jj).airpuff(1)) - refTime;
            else
                lickBeforeAir = tbytDat(jj).Lick(find(tbytDat(jj).Lick > tbytDat(jj).evtOff, 1, 'first'));
                faChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeAir)==1, tbytDat(jj).LickChunk, 'UniformOutput', false));
                faChunkLicks = tbytDat(jj).LickChunk{faChunkI};
                tbytDat(jj).faLicks = faChunkLicks - refTime;
            end

            if ~isempty(tbytDat(jj).airpuff)
                if sum(tbytDat(jj).Lick > tbytDat(jj).airpuff(1))>0
                    % find the first lick after airpuff delivery
                    lickAfterAir = tbytDat(jj).Lick(find(tbytDat(jj).Lick > tbytDat(jj).airpuff(1), 1, 'first'));
                    postAirChunkI = cell2mat(cellfun(@(a) sum(a==lickAfterAir)==1, tbytDat(jj).LickChunk, 'UniformOutput', false));
                    postAirChunkLicks = tbytDat(jj).LickChunk{postAirChunkI};
                    tbytDat(jj).postAirLicks = postAirChunkLicks(postAirChunkLicks>tbytDat(jj).airpuff(1)) - refTime;
                end
            else
                tbytDat(jj).postAirLicks = []; 
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

    function dat = tbytDatRewardPunishI(dat)
        rewardTrI = zeros(1, length(dat));
        punishTrI = zeros(1, length(dat));

        if isfield(dat, 'pos_rwd_tr')
            rewardTrI = rewardTrI+[dat.pos_rwd_tr];
            if isfield(dat, 'pos_oms_tr')
                rewardTrI = rewardTrI+[dat.pos_oms_tr];
                if isfield(dat, 'pos_pns_tr')
                    rewardTrI = rewardTrI+[dat.pos_pns_tr];
                end
            end
        end

        if isfield(dat, 'neg_rwd_tr')
            punishTrI = punishTrI+[dat.neg_rwd_tr];
            if isfield(dat, 'neg_oms_tr')
                punishTrI = punishTrI+[dat.neg_oms_tr];
                if isfield(dat, 'neg_pns_tr')
                    punishTrI = punishTrI+[dat.neg_pns_tr];
                end
            end
        end

        for i = 1:length(dat)
            dat(i).rewardTrI = rewardTrI(i);
            dat(i).punishTrI = punishTrI(i);

        end
    end



end