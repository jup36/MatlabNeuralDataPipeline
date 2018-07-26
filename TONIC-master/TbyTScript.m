%% TRIAL BY TRIAL ANALYSIS SCRIPT

%% Calculate behavioral variables
NumSessions = size(PopData.session,2);
k=1; clear sumBehavior;
validIds = PopData.validBehavSess;

% licking is stored in PopData.session(i).events.EL.ts

for m=1:numel(validIds)

    i = validIds(m);
    
    if numel(PopData.session(i).events.EL.ts)>100
        % instantaneous ILI
        [isi] = TNC_QuantISI(PopData.session(i).events.EL.ts);
        sumBehavior.linTimes        = isi.hist.linTimes;
        sumBehavior.linCount(k,:)   = isi.hist.linCount;
        figure(11); plot(sumBehavior.linCount(k,:),'k'); drawnow;
        k = k+1;
    end

    disp(['Session number: ' num2str(i)]);

    PopData.session(i).behavior.numTrials = size(PopData.session(i).behavior.raster.trial,2);
    
    for j = 1:size(PopData.session(i).behavior.raster.trial,2)

        tmp                 = find(PopData.session(i).behavior.raster.trial(j).ts<2000);
        
        if numel(tmp) > 0
            
            licksToTest     = PopData.session(i).behavior.raster.trial(j).ts(tmp);

            allPostCSLicks  = find(PopData.session(i).behavior.raster.trial(j).ts>0);              
            allAntiCSLicks  = find(licksToTest>0);              
            allPreCSLicks   = find(PopData.session(i).behavior.raster.trial(j).ts<0);

            PopData.session(i).behavior.anticipLicks(j) = numel(allAntiCSLicks) - numel(allPreCSLicks);
            
            % Find the first lick after the cs
            if numel(allPostCSLicks) > 0
                PopData.session(i).behavior.flCSon(j)       = PopData.session(i).behavior.raster.trial(j).ts(allPostCSLicks(1,1));
            else
                PopData.session(i).behavior.flCSon(j)       = 0;
            end
        
            numPL           = numel(allPostCSLicks);

            if numel(allPostCSLicks) > 2

                % For all postLicks try looking for whether they were at a criterion rate (62.63 - 188.65 ms = +/-2sd)
                tmpPostLick = PopData.session(i).behavior.raster.trial(j).ts(allPostCSLicks);
                postLickIsi = tmpPostLick(2:numPL) - tmpPostLick(1:numPL-1);
                critLick    = find(postLickIsi<188.65,1);

                if numel(critLick) == 1
                    PopData.session(i).behavior.fCritLick(j) = PopData.session(i).behavior.raster.trial(j).ts(allPostCSLicks(critLick));
                else
                    PopData.session(i).behavior.fCritLick(j)= 0;                    
                end

            else

                PopData.session(i).behavior.flCSon(j)   = 0;
                PopData.session(i).behavior.fCritLick(j)= 0;

            end
            
        else

                PopData.session(i).behavior.anticipLicks(j) = 0;
                PopData.session(i).behavior.flCSon(j)       = 0;
                PopData.session(i).behavior.fCritLick(j)    = 0;

        end
        
        clear licksToTest allPostCSLicks allAntiCSLicks allPreCSLicks tmpPostLick postLickIsi;

    end
    
                
end

%% Using the newly extracted trial by trial lick counts build up an average behavior during learning and extinction
NumSessions = size(PopData.session,2);
k=1; l=1;
validIds = PopData.validBehavSess;

for m=1:numel(validIds)

    i = validIds(m);

    if strcmp(PopData.session(i).sessClass,'learning')

        for j = 1:numel(PopData.session(i).behavior.anticipLicks)
            allTbyTBehavior.acq.aLbT(j,k)    = PopData.session(i).behavior.anticipLicks(j);            
        end
        k = k+1;
        
    elseif strcmp(PopData.session(i).sessClass,'extinction')

        for j = 1:numel(PopData.session(i).behavior.anticipLicks)
            allTbyTBehavior.ext.aLbT(j,l)    = PopData.session(i).behavior.anticipLicks(j);            
        end
        l = l+1;
        
    else
        
        % this doesn't exist just trying to make the code more clear

    end
    
end

%% Output of the IGOR fit to the peak of the interlick interval
% CurveFit/M=2/W=0 gauss, wave0[75,180]/D
%   Fit converged properly
%   Curve fit with data subrange:
% 	wave0[75,180]
%   fit_wave0= W_coef[0]+W_coef[1]*exp(-((x-W_coef[2])/W_coef[3])^2)
%   W_coef={0.0015453,0.0062659,125.64,31.504}
%   V_chisq= 5.04043e-06;V_npnts= 106;V_numNaNs= 0;V_numINFs= 0;
%   V_startRow= 75;V_endRow= 180;
%   W_sigma={9.86e-05,9.18e-05,0.179,0.628}
%   Coefficient values ± one standard deviation
%   	y0   	=0.0015453 ± 9.86e-05
%   	A    	=0.0062659 ± 9.18e-05
%   	x0   	=125.64 ± 0.179
%   	width	=31.504 ± 0.628


%% define a series of windows and look at spikes within those windows to get trial by trial responses

    currWindow  = respProps.fw(i).X(1,1):1:respProps.fw(i).X(1,2);

    sessNum    =   allCSaligned.sess(i);
    unitNum    =   allCSaligned.unit(i);
    
    % store the current raster in a local variable
    
    % go through all trials and count spikes in the phasic response window


%% pairwise comparisons