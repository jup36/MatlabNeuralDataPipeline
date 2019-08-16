function  trj = interplh(trj,trjlh)
% interpolates a timeseries (trj) based on it's likelihood/a confidence estimate (trjlh)  
x = 1:length(trj);
trj(trjlh<.5)=nan;
%trjwnan = trj; 
%figure; hold on; 
trj(isnan(trj)) = interp1(x(~isnan(trj)),trj(~isnan(trj)),x(isnan(trj)),'spline');
%plot(trj)
%plot(trjwnan)
%hold off; 
end
