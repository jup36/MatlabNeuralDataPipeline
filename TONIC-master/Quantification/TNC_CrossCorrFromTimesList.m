function [TNCxcorr_x,TNCxcorr_y] = TNC_CrossCorrFromTimesList(times1,times2,maxlag)

lastTime = max([max(times1) max(times2)]);

delta1 = zeros(1,lastTime);
delta2 = zeros(1,lastTime);

delta1(round(times1)) = 1;
delta2(round(times2)) = 1;

TNCxcorr_y  = xcorr(delta1,delta2,maxlag,'coeff');
TNCxcorr_x  = -maxlag:maxlag;
