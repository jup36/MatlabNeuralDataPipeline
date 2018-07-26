function [peak] = TNC_QuantPeak(data,threshold,peakVal,win,slope,minSpace)
% FUNCTION DETAILS: input data must be a vector. From that, parameterized by the threshold, window, and sign of slope this function extracts all peaks that it finds.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: dudmanlab.org/projects.html
% _________________________________________________________________________
% 

i=win+1;
k=0; % number of peaks
numSamples = numel(data);

% walk through the data looking for threshold crossings with positive slopes
while i<numSamples-win % give a little space to look for local properties

    if data(i) > threshold
        
        % is the slope the correct direction?
        tmp = polyfit([0:10]',data(i-5:i+5),1);
        localSlope = tmp(1,1);
        
        if slope==1
            if localSlope > 0
                if max(data(i:i+win))>=peakVal
                    truePeak = 1;
                else
                    truePeak = 0;
                end
            else
                truePeak = 0;
            end
        else
            if localSlope < 0
                if min(data(i:i+win))<=peakVal
                    truePeak = 1;
                else
                    truePeak = 0;
                end
            else
                truePeak = 0;
            end            
        end        
        
        if truePeak

            % increment the peak counter
            k = k+1;

            % store this point as the simple threshold detection point
            peak.global.begin(k) = i;
            peak.global.slope(k) = localSlope;

            % look for first place where data exceeded some local run properties
            recentData      = data(i-win:i);
%             avg             = mean(data(i-win:i-3));
%             err             = std(data(i-win:i-3));
%             localThresh     = threshold ./ 2;

%             % look forward to find the first point where it rose above threshold
%             firstUp = find(recentData<localThresh,1,'last');
%             if numel(firstUp)==0
%                 peak.local.begin(k) = i;
%             else
%                 peak.local.begin(k) = i-win+firstUp;
%             end
            
            % look forward to find the last point where it falls back below threshold
            futureData      = data(i+1:numSamples);
            firstDownSimp   = find(futureData<threshold,1,'first');
%             firstDown       = find(futureData<localThresh,1,'first');
            
            peak.global.end(k)  = firstDownSimp+i+1;
%             peak.local.end(k)   = firstDown+i+1;

%             figure(10); plot(data(peak.local.begin(k):peak.local.end(k))); drawnow;
            
            % increment i to continue searching
            i = peak.global.end(k)+minSpace; %enforce a minimum spacing
            
%             disp(['Found peak ' num2str(k) ' between ' num2str(peak.global.begin(k)) ' and ' num2str(peak.global.end(k))]);

        else

            % increment, this is not a rising peak
            i = i+1;

        end

    else
        
        % increment until a threshold is found
        i = i+1;

    end

    truePeak = 0;
    
end

disp(['Found ' num2str(k) ' peaks.']);
peak.numPeaks = k;