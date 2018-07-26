% Compute pairwise cross correlations normalized by shuffled data
function [TNCxcorr_x,TNCxcorr_y,TNCxcorr_y_normed] = TNC_CrossCorrNormed(times1,times2,maxlag,numPerm)

    % first get the raw cross correlation
    lastTime = max([max(times1) max(times2)]);

    delta1 = zeros(1,lastTime);
    delta2 = zeros(1,lastTime);

    delta1(round(times1)) = 1;
    delta2(round(times2)) = 1;

    TNCxcorr_y  = xcorr(delta1,delta2,maxlag,'coeff');
    TNCxcorr_x  = -maxlag:maxlag;

    % get the list of isis
    isi1 = diff(times1);
    isi2 = diff(times2);

    for i=1:numPerm

        times2n(1) = times2(1);
        permTimes = randperm(numel(times2)-1);
        for k=2:numel(permTimes)
            times2n(k) = times2n(k-1)+isi2(permTimes(k));
        end

        lastTime = max([max(times1) max(times2n)]);

        delta1 = zeros(1,lastTime);
        delta2 = zeros(1,lastTime);

        delta1(round(times1)) = 1;
        delta2(round(times2n)) = 1;

        TNCxcorr_Ny(i,:)    = xcorr(delta1,delta2,maxlag,'coeff');
        TNCxcorr_Nx         = -maxlag:maxlag;

    end

    TNCxcorr_y_normed.mn = mean(TNCxcorr_Ny,1);
    TNCxcorr_y_normed.sd = std(TNCxcorr_Ny,[],1);
