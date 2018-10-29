function [ meanmat, stdmat, semmat ] = meanstdsem( inputmat )
%This function receives a single or multiple column(s), and returns mean,
% std, and sem of each column vector.  

numbcol = size(inputmat,2);     % # of columns 
meanmat = nanmean(inputmat,1);      % get mean 
stdmat = nanstd(inputmat,0,1);      % get std

for i = 1:numbcol
    nmat(1,i) = size(find(isnan(inputmat(:,i))==0),1);      % get the valid # of elements
end

semmat = stdmat./sqrt(nmat);

end

