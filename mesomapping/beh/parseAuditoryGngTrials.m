function tbytDat = parseAuditoryGngTrials(tbytDat)

tbytDat = tbytDatRewardPunishI(tbytDat);

chunkCutoff = 0.4; % 0.4 s licks that occur within this cutoff from one another are chunked together
for tt = 1:length(tbytDat)
    if ~isempty(tbytDat(tt).Lick)
        refTime = tbytDat(tt).evtOn;

        % Chunk licks
        tbytDat(tt).LickChunk = chunkTimestamps(tbytDat(tt).Lick, chunkCutoff);

        % ITI (find lick chunks where every element belongs to ITI)
        itiChunkI = cell2mat(cellfun(@(a) sum(a < tbytDat(tt).evtOn)==length(a), tbytDat(tt).LickChunk, 'UniformOutput', false));
        tbytDat(tt).itiLickChunk = cellfun(@(a) a-refTime, tbytDat(tt).LickChunk(itiChunkI), 'UniformOutput', false);

        %% classify types of licks in go/no-go trials
        % Hit
        if (tbytDat(tt).pos_rwd_tr == 1 && ~isempty(tbytDat(tt).water)) || ...
                (tbytDat(tt).pos_oms_tr == 1 && sum(tbytDat(tt).Lick>tbytDat(tt).evtOff)>0) || ...
                (tbytDat(tt).pos_pns_tr == 1 && sum(tbytDat(tt).Lick>tbytDat(tt).evtOff)>0)
            % find the last lick before water delivery that likely triggered the reward
            if ~isempty(tbytDat(tt).water)
                lickBeforeWater = tbytDat(tt).Lick(find(tbytDat(tt).Lick <= tbytDat(tt).water(1), 1, 'last'));
                if isempty(lickBeforeWater)
                    lickBeforeWater = tbytDat(tt).Lick(1); 
                    hitChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeWater)==1, tbytDat(tt).LickChunk, 'UniformOutput', false));
                    hitChunkLicks = tbytDat(tt).LickChunk{hitChunkI};
                    tbytDat(tt).hitLicks = hitChunkLicks(1) - refTime;
                else
                    hitChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeWater)==1, tbytDat(tt).LickChunk, 'UniformOutput', false));
                    hitChunkLicks = tbytDat(tt).LickChunk{hitChunkI};
                    tbytDat(tt).hitLicks = hitChunkLicks(hitChunkLicks<=tbytDat(tt).water(1)) - refTime;
                end
            else
                lickBeforeWater = tbytDat(tt).Lick(find(tbytDat(tt).Lick > tbytDat(tt).evtOff, 1, 'first'));
                hitChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeWater)==1, tbytDat(tt).LickChunk, 'UniformOutput', false));
                hitChunkLicks = tbytDat(tt).LickChunk{hitChunkI};
                tbytDat(tt).hitLicks = hitChunkLicks - refTime;
            end

            % find licks after the stim offset
            lickAfterevtOffI = cell2mat(cellfun(@(a) sum(a>=tbytDat(tt).evtOff)==length(a), tbytDat(tt).LickChunk, 'UniformOutput', false));
            tbytDat(tt).postStimChunk = cellfun(@(a) a-refTime, tbytDat(tt).LickChunk(lickAfterevtOffI), 'UniformOutput', false);

            if ~isempty(tbytDat(tt).water)
                % Consumptive licks
                if sum(tbytDat(tt).Lick > tbytDat(tt).water(1))>0
                    % find the first lick after water delivery
                    lickAfterWater = tbytDat(tt).Lick(find(tbytDat(tt).Lick > tbytDat(tt).water(1), 1, 'first'));
                    consumeChunkI = cell2mat(cellfun(@(a) sum(a==lickAfterWater)==1, tbytDat(tt).LickChunk, 'UniformOutput', false));
                    consumeChunkLicks = tbytDat(tt).LickChunk{consumeChunkI};
                    tbytDat(tt).consumeLicks = consumeChunkLicks(consumeChunkLicks>tbytDat(tt).water(1)) - refTime;
                end
            else
                tbytDat(tt).consumeLicks = []; 
            end
        % False Alarm
        elseif (tbytDat(tt).neg_rwd_tr == 1 && ~isempty(tbytDat(tt).water)) || ...
                (tbytDat(tt).neg_oms_tr == 1 && sum(tbytDat(tt).Lick>tbytDat(tt).evtOff)>0) || ...
                (tbytDat(tt).neg_pns_tr == 1 && sum(tbytDat(tt).Lick>tbytDat(tt).evtOff)>0)
            % find the last lick before airpuff that likely triggered the
            % airpuff
            if ~isempty(tbytDat(tt).airpuff)
                if tbytDat(tt).Lick(1) <= tbytDat(tt).airpuff(1)
                    lickBeforeAir = tbytDat(tt).Lick(find(tbytDat(tt).Lick < tbytDat(tt).airpuff(1), 1, 'last'));
                    faChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeAir)==1, tbytDat(tt).LickChunk, 'UniformOutput', false));
                    faChunkLicks = tbytDat(tt).LickChunk{faChunkI};
                    tbytDat(tt).faLicks = faChunkLicks(faChunkLicks<=tbytDat(tt).airpuff(1)) - refTime;
                end
            else
                lickBeforeAir = tbytDat(tt).Lick(find(tbytDat(tt).Lick > tbytDat(tt).evtOff, 1, 'first'));
                faChunkI = cell2mat(cellfun(@(a) sum(a==lickBeforeAir)==1, tbytDat(tt).LickChunk, 'UniformOutput', false));
                faChunkLicks = tbytDat(tt).LickChunk{faChunkI};
                tbytDat(tt).faLicks = faChunkLicks - refTime;
            end

            if ~isempty(tbytDat(tt).airpuff)
                if sum(tbytDat(tt).Lick > tbytDat(tt).airpuff(1))>0
                    % find the first lick after airpuff delivery
                    lickAfterAir = tbytDat(tt).Lick(find(tbytDat(tt).Lick > tbytDat(tt).airpuff(1), 1, 'first'));
                    postAirChunkI = cell2mat(cellfun(@(a) sum(a==lickAfterAir)==1, tbytDat(tt).LickChunk, 'UniformOutput', false));
                    postAirChunkLicks = tbytDat(tt).LickChunk{postAirChunkI};
                    tbytDat(tt).postAirLicks = postAirChunkLicks(postAirChunkLicks>tbytDat(tt).airpuff(1)) - refTime;
                end
            else
                tbytDat(tt).postAirLicks = []; 
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