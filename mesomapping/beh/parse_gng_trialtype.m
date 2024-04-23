function tt = parse_gng_trialtype(stimopts)
%This function takes stimopts from a gonogo session, and parse the trial
% types. 

%% parse positive stim 
tt.pos_rwd_tr = zeros(length(stimopts.rewarded_stim), 1);
tt.pos_oms_tr = zeros(length(stimopts.rewarded_stim), 1);
tt.pos_pns_tr = zeros(length(stimopts.rewarded_stim), 1);

posTtLogic = stimopts.proba_positive_stim > 0; 
posTr_outcome = [find(stimopts.stim_type==1), stimopts.outcome_positive_stim];

if posTtLogic(1) % reward 
    posTr_rwd_outcomeI = posTr_outcome(:, 2) == 1;
    tt.pos_rwd_tr(posTr_outcome(posTr_rwd_outcomeI, 1)) = 1;
end

if posTtLogic(2) % omission
    posTr_oms_outcomeI = posTr_outcome(:, 2) == 0;
    tt.pos_oms_tr(posTr_outcome(posTr_oms_outcomeI, 1)) = 1;
end

if posTtLogic(3) % punishment
    posTr_pns_outcomeI = posTr_outcome(:, 2) == -1;
    tt.pos_pns_tr(posTr_outcome(posTr_pns_outcomeI, 1)) = 1;
end

%% parse negative stim 
tt.neg_rwd_tr = zeros(length(stimopts.punished_stim), 1);
tt.neg_oms_tr = zeros(length(stimopts.punished_stim), 1);
tt.neg_pns_tr = zeros(length(stimopts.punished_stim), 1);

negTtLogic = stimopts.proba_negative_stim > 0; 
negTr_outcome = [find(stimopts.stim_type==2), stimopts.outcome_negative_stim];

if negTtLogic(1) % reward 
    negTr_rwd_outcomeI = negTr_outcome(:, 2) == 1;
    tt.neg_rwd_tr(negTr_outcome(negTr_rwd_outcomeI, 1)) = 1;
end

if negTtLogic(2) % omission
    negTr_oms_outcomeI = negTr_outcome(:, 2) == 0;
    tt.neg_oms_tr(negTr_outcome(negTr_oms_outcomeI, 1)) = 1;
end

if negTtLogic(3) % punishment
    negTr_pns_outcomeI = negTr_outcome(:, 2) == -1;
    tt.neg_pns_tr(negTr_outcome(negTr_pns_outcomeI, 1)) = 1;
end

end