function  [intTrj,fstPt] = interplh(trj,trjlh,lm, dealFirstNaNs)
% interpolates a timeseries (trj) based on it's likelihood (trjlh), enforce the Trjs to be within the frame dimension (lm)  
%trj = tmpTrj.fg1Y;

%trj = tmpTrj.js1X; trj1h = tmpTrj.js1lh; lm = wd; dealFirstNaNs = true;
intTrj = nan(length(trj),1); 
x = 1:length(trj);

trj(trjlh<.9)=nan;
fstPt = NaN;

%trjSpl = trj; 
%trjSpl(trjlh<.1)=nan;

%trjPch = trj; 
%trjPch(trjlh<.1)=nan;

%sum(isnan(trj))/length(trj); 

if sum(isnan(trj(end-4:end)))==0 && sum(isnan(trj))/length(trj)<.4 % interpolation would be problematic with NaNs at the end
    
    if dealFirstNaNs % initial NaNs cannot be interpolated (e.g. joystick before/during positioning) - just interpolate subsequent to the 1st non-NaN value
        fstPt = find(trjlh>.9,1,'first'); 
        tempTrj = trj(fstPt:end); 
        if sum(isnan(tempTrj))>0
            x = 1:length(tempTrj); 
            tempTrj(isnan(tempTrj)) = interp1(x(~isnan(tempTrj)),tempTrj(~isnan(tempTrj)),x(isnan(tempTrj)),'pchip');
        end
        
        trj(fstPt:end) = tempTrj; 
        trj = min(trj,lm); 
        trj = max(trj,0); 
        trj(1:fstPt) = 0;
        intTrj = trj; 
    %trjwnan = trj;
    %figure; hold on;
    %trj(isnan(trj)) = interp1(x(~isnan(trj)),trj(~isnan(trj)),x(isnan(trj)),'spline');
    %trjSpl(isnan(trj)) = interp1(x(~isnan(trj)),trj(~isnan(trj)),x(isnan(trj)),'spline');
    %trjPch(isnan(trj)) = interp1(x(~isnan(trj)),trj(~isnan(trj)),x(isnan(trj)),'pchip');
    else % initial NaNs are not really an issue for body parts
        trj(isnan(trj)) = interp1(x(~isnan(trj)),trj(~isnan(trj)),x(isnan(trj)),'pchip');
        if sum(trj>lm)+sum(trj<0)==0 % enforce the trajectories within the frame
            intTrj = trj;
        end    
        fstPt = NaN;
    %hold on; 
    %plot(trj); plot(trjSpl); plot(trjPch); 
    %plot(trj)
    %plot(trjwnan)
    %hold off;
    end
    

end

end