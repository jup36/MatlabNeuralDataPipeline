function [ afterfiltermat ] = twostdfilter( beforefiltermat )
%This fuction is to remove outlier elements 'OF EACH COLUMN' that are out of the +- 2
%   standard deviation range. The input 'beforefiltermat' can be either
%   single or multiple column matrix. 

meanmat = nanmean(beforefiltermat,1);                                       % get the mean 
stdmat = nanstd(beforefiltermat,0,1);                                       % get the standard deviation 
upperboundary = meanmat + 2.*stdmat;                                        % set the upper boundary
lowerboundary = meanmat - 2.*stdmat;                                        % set the lower boundary
afterfiltermat = nan(size(beforefiltermat,1), size(beforefiltermat,2));     % preallocate the output matrix filled with NaN

for i = 1:size(beforefiltermat,2)
    validelements = find((beforefiltermat(:,i) <= upperboundary(1,i) & (beforefiltermat(:,i) >= lowerboundary(1,i))));
    afterfiltermat(validelements,i) = beforefiltermat(validelements,i);
end

