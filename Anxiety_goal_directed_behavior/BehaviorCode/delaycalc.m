function [ delay ] = delaycalc( onset, offset, timeconstant )
% This function is used to get the delay from each onset to offset

if size(onset,1)==size(offset,1) 
   delay = offset - onset;
else
   error('myApp:argChk', 'The number of onset and offset elements must match!') 
end

delay = delay.*timeconstant;

return

