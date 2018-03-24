% This function is to perform z-score normalization of the raw power
% spectral density values to the mean baseline power spectral
% densities. The input 'rawpsdcell' contains collected baseline psd values
% of all trials (not classified based on the types of trials).
% This function returns z-score normalized psds in the 'z-cell',
% but returns NaN if the raw psds were NaNs or there was only one trial in
% the raw psd cell which leads to standard deviation = 0 (to avoid dividing with zero).
% This function also returns the mean and std of the collected baseline psds.
% Modified on 03/12/13 to deal with 

function [ z_cell, basepsdmean, basepsdstd ] = psdnormalization(rawpsdcell, validbasecell, t, f)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

basemat = cell2mat(validbasecell);

if isempty(basemat) == 0        % this means that there's at least one valid trial. 
    
    basepsdmean = nanmean(basemat,1);
    basepsdstd = nanstd(basemat,0,1);
    
    z_cell = cell(size(rawpsdcell,1),size(rawpsdcell,2));
    
    if any(basepsdstd == 0)             % if any element in the std array is zero this function must return NaN to avoid dividing with zero (infinite)
        for i = 1:size(rawpsdcell,1)
            for ii = 1:size(rawpsdcell,2)
                z_cell{i,ii} = NaN;
            end
        end
    else
        for i = 1:size(rawpsdcell,1)
            for ii = 1:size(rawpsdcell,2)
                tempnanidx = isnan(rawpsdcell{i,ii});
                if sum(tempnanidx(:)) == 0
                    z_cell{i,ii} = zeros(size(rawpsdcell{i,ii},1),size(rawpsdcell{i,ii},2));
                    for j = 1:size(rawpsdcell{i,ii},1)         % # of rows in each element of the cell, must be 20
                        z_cell{i,ii}(j,:) = (rawpsdcell{i,ii}(j,:)-basepsdmean)./basepsdstd;
                    end
                else % if the rawpsdcell{i,ii} is NaN;
                    z_cell{i,ii} = NaN;
                end
            end
        end
    end
    
else                            % this means that there's no valid trial at all.
    basepsdmean = NaN(1,size(f,2));         % put NaNs for the empty matrix. 
    basepsdstd = NaN(1,size(f,2));          % put NaNs for the empty matrix. 
    for i = 1:size(rawpsdcell,1)
        for ii = 1:size(rawpsdcell,2)
            z_cell{i,ii} = NaN;
        end
    end
end