function  [intTrj,fstPt] = interplh(trj,trjlh,lm, zeroFirstNaNs)
% interpolates a timeseries (trj) based on it's likelihood (trjlh), enforce the Trjs to be within the frame dimension (lm)
%trj = tmpTrj.fg1Y;

%trj = tmpTrj.js1X; trjlh = tmpTrj.js1lh; lm = wd; zeroFirstNaNs = true;
intTrj = nan(length(trj),1);
x = 1:length(trj);

trj(trjlh<.9)=nan;
fstPt = NaN;

if zeroFirstNaNs
    isnanThres = .9; % in case for Joystick, there's a lot of NaNs expected while Js positioning
else
    isnanThres = .4; 
end

if sum(isnan(trj))/length(trj)<isnanThres % interpolation would be problematic with NaNs at the end
    valPts = strfind(num2str(trjlh>.9)','111');
    if sum(isnan(trj(end-4:end)))==0 % deal with last NaNs
        lastValPt = length(trj);
    elseif ~isempty(valPts)
        lastValPt = valPts(end)+2; % last valid point
    end
   
    fstPt = find(trjlh>.9,1,'first');
    tempTrj = trj(fstPt:lastValPt); % the valid portion of the trajectory
    tempTrjIdx = 1:length(trj)>=fstPt & 1:length(trj)<=lastValPt; % index for valid points
    if sum(isnan(tempTrj))>0
        % intrapolate points between the first and last valid points
        tempTrj(isnan(tempTrj)) = interp1(x(~isnan(tempTrj)),tempTrj(~isnan(tempTrj)),x(isnan(tempTrj)),'pchip');
    end
    trj(tempTrjIdx) = tempTrj;
    % extrapolate if each of the invalid portions at both sides is less than .1 of the full trj length
    if fstPt < length(trj)*.1 && lastValPt > length(trj)*.9
        trj(tempTrjIdx==0) = interp1(x(tempTrjIdx),trj(tempTrjIdx),x(tempTrjIdx==0),'pchip','extrap');
    end
    % put zeros for first NaNs, which often is the case for joystick trajectories
    if zeroFirstNaNs
        trj(1:fstPt) = 0;
    end
    trj = min(trj,lm);
    trj = max(trj,0);
    intTrj = trj;
end

end

