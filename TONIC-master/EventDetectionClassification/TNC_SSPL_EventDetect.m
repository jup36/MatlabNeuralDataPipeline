% FUNCTION DETAILS: This function goes through a single channel of filtered recording data and looks for threshold crossings. A second stage then tests these threshold crossings according to a template matching heuristic to try to classify significant events.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: josh@dudmanlab.org
% CONTRIBUTIONS: www.dudmanlab.html
%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%
% _________________________________________________________________________
% OUTPUT STRUCTURE ELEMENTS
% events.wfs
% events.inds
% events.sampleRate
% events.channel
% events.snrThresh
% events.winL
% events.winR
% events.numEvents

function [events] = TNC_SSPL_EventDetect(bandPassedData,sampleRate,snrThresh,timeLockOut)

 
%% EVENT DETECTION
tstart = tic;
 
    peakTime = 1;
 
    events.sampleRate = sampleRate;
    samples = numel(bandPassedData);
    
% find the rms values of the data to set the SNR threshold
    % standard estimator
    stdRawData = std(bandPassedData);
    
    % robust estimator
    avgRawData = mean(median(abs(bandPassedData)));
    stdRawData = median(abs(bandPassedData-avgRawData)) ./ 0.6745;
    

    MAD = median(abs(bandPassedData - median(bandPassedData)));
    stdRawData1 = MAD ./ 0.6745;
    stdRawData2 = (quantile(bandPassedData, 0.75) -  quantile(bandPassedData, 0.25))/1.35;
%     disp('  ');
%     disp(['max(bandPassedData)=' num2str(max(bandPassedData)) ' min(bandPassedData)=' num2str(min(bandPassedData))]);
    N = round((max(bandPassedData) - min(bandPassedData))/(stdRawData1/20));
%     disp(['N=' num2str(N)]);
    disp(['avgRawData=' num2str(avgRawData) ' stdRawData=' num2str(stdRawData) ' stdRawData1=' num2str(stdRawData1) ' stdRawData2=' num2str(stdRawData2)]);
    spikeFreeData = bandPassedData(abs(bandPassedData) < abs(snrThresh)*stdRawData1);
    sigma_fit = get_spike_free_sample_stddev(spikeFreeData, N);

    if snrThresh < 0
        threshold = -snrThresh;
    else
        threshold = stdRawData1.*snrThresh;
%       threshold = avgRawData.*snrThresh;
    end
    events.ampThresh = threshold;
%     disp(['threshold=' num2str(threshold) ' sigma_fit=' num2str(sigma_fit) ' avgRawData=' num2str(avgRawData) ' stdRawData=' num2str(stdRawData)]);
%     disp(' ');
%     disp(' ');
%     disp(' ');

%   disp(['sigma=' num2str(sigma)]);

%     disp(['SD calculated to be: ' num2str(stdRawData) ' ... giving a threshold of: ' num2str(threshold)]);
 
% display the threshold choice
%     figure(2); clf;
    xValsForThresh = 1:1000:length(bandPassedData);
    yValsForThresh = ones(1,length(xValsForThresh)) .* threshold;
%     plot(1:1:length(bandPassedData),bandPassedData(1,:),'b',xValsForThresh,yValsForThresh,'r-',xValsForThresh,-yValsForThresh,'r-')
%     drawnow;
    
% first stage is to run find on the array of data values
    if threshold>0
        allSupraThreshIndices = find(bandPassedData<-threshold);
    else
        allSupraThreshIndices = find(bandPassedData<threshold);        
    end

%     disp(sprintf('First pass detected events: %g',length(allSupraThreshIndices)));
%     plot(1:1:length(bandPassedData),bandPassedData(1,:),'b',allSupraThreshIndices,bandPassedData(allSupraThreshIndices),'r.');
%     drawnow;
    
    if numel(allSupraThreshIndices)==0
 
        candidateEvents = [];
        events.inds = [];
        events.ts   = [];
        
    else
 
    % eliminate any crossings without a return below threshold
        candidateEvents(1,1) = allSupraThreshIndices(1,1);
        j=2;
        for i=2:length(allSupraThreshIndices)
            if  allSupraThreshIndices(1,i)<(samples-60) && allSupraThreshIndices(1,i)>(60)
                if (allSupraThreshIndices(1,i) ~= allSupraThreshIndices(1,i-1) + 1) & (allSupraThreshIndices(1,i) > allSupraThreshIndices(1,i-1) + (timeLockOut * 30) )
                    if peakTime==1
                        % should reach a minima by 0.4 ms after threshold cross
                        pkLoc = find(bandPassedData(1,allSupraThreshIndices(1,i):allSupraThreshIndices(1,i)+12) == min(bandPassedData(1,allSupraThreshIndices(1,i):allSupraThreshIndices(1,i)+12)),1);
                        % should also change by at least 2*Threshold within +/- 0.27ms (8 samples) 
                        pkMin = bandPassedData(1,allSupraThreshIndices(1,i)+pkLoc);
%                         twoT = pkMin+(2.*threshold);
%                         dmpCriterion = find(bandPassedData(1,allSupraThreshIndices(1,i)+(pkLoc-1):allSupraThreshIndices(1,i)+(pkLoc+8))>=twoT);
                        dmpCriterion = find(bandPassedData(1,allSupraThreshIndices(1,i)+(pkLoc-12):allSupraThreshIndices(1,i)+(pkLoc+12))>=(0));
                        if numel(pkLoc)==1 && numel(dmpCriterion)>0
                            candidateEvents(1,j) = allSupraThreshIndices(1,i)+pkLoc-1; 
                            j=j+1;
                        else
%                             disp('not a real spike');
                        end
                    else
                        candidateEvents(1,j) = allSupraThreshIndices(1,i)-2;
                        j=j+1;
                    end
 
                end
            end  
        end
 
        events.inds = candidateEvents;
        events.ts   = candidateEvents./sampleRate;
 
    end
    
telapsed = toc(tstart);
 
    numberOfEvents = length(candidateEvents);    
 
    disp(sprintf('     EVENT DETECT | Number of Events: %g | Elapsed time (seconds): %g',numberOfEvents,telapsed));
 
%-------------------------------------------------------------------------------

function sigma = get_spike_free_sample_stddev(data, N)
%     disp(['median(data)=' num2str(median(data))]);
%   data        = data - median(data);
%   x_thr       = min(max(data), -min(data));
%   data        = data(data <= x_thr & data >= -x_thr);
    std_est     = cap_fit(data, N);
    scale_data  = max(data)-min(data);
    F           = normcdf(data, 0., std_est);
    X0_9987_max = max(data(F <= 0.9987));
    X0_9987_min = min(data(F >= 0.9987));
    alpha       = 3.1564;
%   alpha       = 3.;
%     disp(['numel(data)=' num2str(numel(data)) ' X0_9987_max=' num2str(X0_9987_max) ' X0_9987_min=' num2str(X0_9987_min) ' std_est*alpha=' num2str(std_est*alpha)]);
    data1       = data;
    iter = 0;
    while min(X0_9987_min,X0_9987_max) > std_est*alpha
        iter = iter + 1;
        disp(['    Iter=' num2str(iter)]);
        x_thr   = 0.99*min(max(data1), -min(data));
        data1   = data1(data1 <= x_thr & data1 >= -x_thr);
        std_est = cap_fit(data1, N);
        F       = normcdf(data1, 0., std_est);
        X0_9987_min = min(data1(F >=0.9987));
        X0_9987_max = max(data(F <= 0.9987));
        disp(['    numel(data)=' num2str(numel(data1)) ' X0_9987=' num2str(min(X0_9987_min+X0_9987_max)) ' std_est*alpha=' num2str(std_est*alpha) ]);
    end
    sigma       = std_est*scale_data;

%-------------------------------------------------------------------------------

% Estimate the standard deviation of the noise sample
%
% See: Pramodsingh H.Thakur et al.
%      Automated optimal detection and classification of neural
%      action potentials in extra-cellular recordings
%      Journal of Neuroscience Methods 162 (2007) 364.376
function sigma = cap_fit(data0, N)
    scale_data = max(data0)-min(data0);
    data = data0/scale_data;
    [y0, x]=hist(data, N);
    n0 = find(y0 == max(y0)); % coord of maximum
    nbeg = n0-100;
    nend = n0+100;
    inds = [(n0-60):(n0+60)];
 %  dn = min(nend - n0, n0 - nbeg);
 %  nbeg = n0 - dn;
 %  nend = n0 + dn;
    y0 = y0(inds);
    x  = x( inds);     
%     disp(['y0=' num2str(y0)]);
%     disp(['x =' num2str(x )]);
    N = numel(inds);               

    MAD = median(abs(data - median(data)));
    sigma_rob = MAD ./ 0.6745;
    sigma_reg = std(data);
    A1 = -1/(2.*sigma_rob^2);
    B1 = - log(sqrt(2.*pi*sigma_rob^2));
    A2 = -1/(2.*sigma_reg^2);
    B2 = - log(sqrt(2.*pi*sigma_reg^2));
%     B1
    zs = A1*x.^2 + B1;
    Integral = sum(y0).*(x(2)-x(1));
%     disp(['N_final=' num2str(N) ' Integral=' num2str(Integral) ' x(1)=' num2str(x(1)) ' x(N)=' num2str(x(N))]);
%     disp(['sigma_robust=' num2str(sigma_rob*scale_data) ' A_robust=' num2str(A1) ' B_robust=' num2str(B1)]);
%     disp(['sigma_reg=   ' num2str(sigma_reg*scale_data) ' A_reg='    num2str(A2) ' B_reg='    num2str(B2)]);

%   plot(x, y0/Integral, 'Color',[1 0 0]); hold on;  % bell-shaped curve
%   plot(x, y, 'Color',[1 0 0]); hold on;  
%   plot(x.^2, y, 'Color',[1 0 0]); hold on;

    % Analysis according to the paper
    y = log(y0/Integral);
    A_paper = (mean(y.*(x.^2)) - mean(y)*mean(x.^2))/(mean(x.^4) - mean(x.^2)*mean(x.^2));
    B_paper =  mean(y) - A_paper*mean(x.^2);
%   disp(['x=' num2str(x)]);
%   disp(['y=' num2str(y)]);
    sigma1_paper = sqrt(-1/2/(A_paper));
    sigma2_paper = exp(-B_paper)/sqrt(2*pi);
    zp = A_paper*x.^2 + B_paper;
%     disp(['sigma1_paper=  ' num2str(sigma1_paper*scale_data) ' sigma2_paper=  ' num2str(sigma2_paper*scale_data) ' A_paper=' num2str(A_paper) ' B_paper=' num2str(B_paper)]);

    % Analysis using polyfit
    p=polyfit(x.^2,y,1);
    A0 = p(1);
    B0 = p(2);
    z0 = A0*x.^2 + B0;
    sigma1_polyfit = sqrt(-1/2/(A0));
    sigma2_polyfit = exp(-B0)/sqrt(2*pi);
%     disp(['sigma1_polyfit=' num2str(sigma1_polyfit*scale_data) ' sigma2_polyfit=' num2str(sigma2_polyfit*scale_data) ' A_polyfit=' num2str(A0) ' B_polyfit=' num2str(B0)]);

    % Producing plot
    PRODUCE_PLOT = 0;
    if PRODUCE_PLOT
        figure(1); clf;
        subplot(10,1,1:8);
        scatter(x.^2, y                  ); hold on;
        plot(   x.^2, z0, 'Color',[1 0 0]); hold on;
%       plot(   x.^2, zp, 'Color',[0 1 0]); hold on;
%       z2 = A2*x.^2 + B2;
%       plot(x.^2, zs, 'Color',[0 .5 .5]); hold on;
%       plotregression(y,x,'Regression');
%       plot(x.^2, z, 'Color',[1 0 0]);
        zoom xon;
        drawnow;
        waitfor(figure(1));
%         disp('ok');
    end

    sigma  = max(sigma1_polyfit, sigma2_polyfit);
 

