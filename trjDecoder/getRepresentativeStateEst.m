function [representState,representStateCut,representIntState] = getRepresentativeStateEst(refState,estState,valTrId)
%refState = stateD;
%estState = estStateCtxMean;
%valTrId = valTrI;
nKv = mode(cell2mat(cellfun(@(a) size(a,1),refState(valTrId),'un',0)));
for rr = 1:size(refState,1)
    for cl = 1:size(refState,2)
        if valTrId(rr,cl)
            distToRef = squeeze(cell2mat(cellfun(@(a) sum(sqrt((a-refState{rr,cl}).^2)),estState(rr,cl,:),'un',0))); % average across ctx resampled trials
            representState{rr,cl} = estState{rr,cl,find(distToRef==min(distToRef),1,'first')};
            tmpState = representState{rr,cl};
            tmpCutState = nan(nKv,50);
            if size(tmpState,2)>=5
                tmpCutState(:,1:min(size(tmpState,2),50)) = tmpState(:,1:min(size(tmpState,2),50));
                representStateCut{rr,cl} = tmpCutState;
                representIntState{rr,cl} = intm(tmpState,50); % original estimated trajectory with interpolation
            end
        end
    end
end
clearvars r c
end