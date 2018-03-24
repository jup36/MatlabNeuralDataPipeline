% This function 'eventpsd' is embedded in the LFP_set_shift_mtspecgramc_eventpsd.
% 'eventpsd' requires the voltage signal (d), timestamp (tstp._),
% which is gained by the script 'LFP_set_shift_mtspecgramc_timestamps', movingwin,
% sv, and parameters for 'mtspecgramc' (params).
% 'eventpsd' returns both raw and normalized power
% spectra values in the structure of psd.
% 'eventpsd' returns NaNs if either of baseline or peri-event voltage
% signal violates the set criterion.

% modified on 10/22/12


function [ current_psdmat, current_rawmat, current_norm_psdmat, current_basemat, current_rawbasemat, tempt, tempf ] = eventpsdAGB( currentd, currenttstp, currentcue, movingwin, params, sv, upperlimit, lowerlimit )
% This function is to perform FFT using mtspecgramc to get power spectral
% density around each meanningful event
% Input: currentd(voltage signal), currenttstp(timestamps), movingwin & params(for mtspecgramc), sv(other parameters)
% Output:
% current_psdmat: raw psds without normalization.
% current_rawmat: raw voltage signal.
% current_norm_psdmat: z-score psds normalized to its own baseline.
% current_basemat: raw psds of the baseline period averaged across the
% time bins.
% current_rawbasemat: raw voltage signal of the baseline period.

% This part has been modified on (03/11/13) to use individualized limits
% that will be fed as inputs of this function.
% meand = mean(currentd,2);
% stdd = std(currentd,0,2);
% upperlimit = meand+2*stdd;      % this has been changed from 4*stdd
% lowerlimit = meand-2*stdd;

if isnan(currenttstp) == 0      % in case the current tstp is valid (non-NaN)
    
    tempd = createdatamatc(currentd,currenttstp,1000,[2.25 2.25]);      % get the temporary data using the createdatamatc, +- .25 are to give extra spacing due to the moving window (if the size of moving windows is 500 ms, start at -250 ms becuase the mtspecgrams takes 500 ms from the starting point).
    tempcue = currentcue(max(find(currentcue <= currenttstp)));               % find the ITI of each trial
    tempbased = createdatamatc(currentd,tempcue+sv.baseendpoint,1000,[2 0]);  % get the temporary base data using the createdatamatc, sv.baseendpoint is set to 9.25 sec considering the moving window
    % plot(tempd); axis([0 length(tempd) -1 1]);
    % plot(tempbased); axis([0 length(tempbased) -1 1]);
    [tempS,tempt,tempf] = mtspecgramc(tempd,movingwin,params);      % run fast fourier transform using mtspecgramc, plot_matrix(tempS,tempt,tempf) can be used for plotting;
    [tempbaseS,tempbaset,tempbasef] = mtspecgramc(tempbased,movingwin,params);
    tempS = tempS(1:size(tempS,1)-1,:);             % remove the last row from the matrix, which was not moving-windowed
    tempbaseS = tempbaseS(1:size(tempbaseS,1)-1,:);             % remove the last row from the matrix, which was not moving-windowed (the final # of row should be 20)
    current_psdmat = zeros(length(tempt)-1,length(tempf));      % preallocate the matrix to contain psd values
    current_norm_psdmat = zeros(length(tempt)-1,length(tempf));
    
    %% this loop is exclude of noisy bins from calculating the mean baseline psd, and this is the old version, in which parts of the bins in the baseline period have been used to calculate the mean and std of the baseline period.
    % bintempbased = zeros(size(tempbaseS,1),sv.timewindow);
    % idxbintempbased = zeros(size(tempbaseS,1),1);                   % use this as the index array put 0 if there's no crossing the threshold
    % for jj = 0:size(tempbaseS,1)-1                                  % # of time bins
    %    if jj == 0
    %        bintempbased(jj+1,:) = tempbased(1,1:sv.timewindow);
    %        idxbintempbased(jj+1,1) = sum((bintempbased(jj+1,:)>upperlimit)|(bintempbased(jj+1,:)<lowerlimit));
    %    else
    %        bintempbased(jj+1,:) = tempbased(1,jj*sv.step+1:(jj*sv.step+1)+sv.timewindow-1);
    %        idxbintempbased(jj+1,1) = sum((bintempbased(jj+1,:)>upperlimit)|(bintempbased(jj+1,:)<lowerlimit));
    %    end
    % end
    
    % remove the last row from the matrix, which was not moving-windowed.
    % current_psdmat = current_psdmat(1:size(current_psdmat,1)-1,:);
    % current_norm_psdmat = current_norm_psdmat(1:size(current_norm_psdmat,1)-1,:);
    % current_basemat = tempmeanbaseS;
    
    %% this is the new criterion for baseline, in which the entire baseline will be thrown out if there's any violation of the upper or lower limit.
    if sum((tempbased>upperlimit)|(tempbased<lowerlimit)) >= 1;     % base is noisy, then remove all
        current_psdmat = NaN;
        current_rawmat = NaN;
        current_norm_psdmat = NaN;
        current_basemat = NaN;
        current_rawbasemat = NaN;
    elseif sum((tempd>upperlimit)|(tempd<lowerlimit)) >= 1;         % signal is noisy, then remove all
        current_psdmat = NaN;
        current_rawmat = NaN;
        current_norm_psdmat = NaN;
        current_basemat = NaN;
        current_rawbasemat = NaN;
    else                                                            % noise free, return base and signal
        tempmeanbaseS = mean(tempbaseS,1);
        tempstdbaseS = std(tempbaseS,0,1);
        current_psdmat = tempS;
        current_rawmat = tempd;
        current_basemat = tempmeanbaseS;
        current_rawbasemat = tempbased;
        %% this loop is to perform normalization to the mean baseline psd
        for u = 1:size(tempS,1)
            current_norm_psdmat(u,:) = (tempS(u,:)-tempmeanbaseS)./tempstdbaseS;
        end
    end
    
else        % in case the current tstp is valid (non-NaN)
    current_psdmat = NaN;
    current_rawmat = NaN;
    current_norm_psdmat = NaN;
    current_basemat = NaN;
    current_rawbasemat = NaN;
    tempt = NaN;
    tempf = NaN;
    
end


end



