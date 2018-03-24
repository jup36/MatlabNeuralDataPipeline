function p = bins2psth(x)
% function p = bins2psth(x)
%
% computes the peri-stimulus time histogram of the binned spike
% train x, where the kth trial is in the kth column of x
%
% p = mean(x,2) is a column vector
%
% p(k) = average firing rate of x at time k

p = sum(x,2)./size(x,2);