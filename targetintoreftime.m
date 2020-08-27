function [targinref] = targetintoreftime(A,timeA,timeRef, varargin)
%This function takes a times series 'A' across a time window ('timeA') and
% plug it into the new time window 'timeRef'. 
% 1-ms time resolution is assumed. 
p = targetintoreftime( A,timeA,timeRef, varargin ); 

Aval = A(timeRef(1)<=timeA&timeA<=timeRef(end)); % valid portion of A to be placed in timeRef

if p.Results.fillin == 0 
    targinref = zeros(size(A,1),size(timeRef,2)); 
elseif isnan(p.Results.fillin)
    targinref = nan(size(A,1),size(timeRef,2)); 
end

tRef1 = find(timeRef>=timeA(1),1,'first'); % 1st point
if ~isempty(tRef1) && ~isempty(Aval)
    targinref(tRef1:tRef1+min(length(targinref),length(Aval))-1) = Aval;
end

function p = targetintoreftime( A,timeA,timeRef, vargs )
% parse input, and extract name-value pairs
default_fillin = 0; % by default do not re-read the raw bin file, if done already

p = inputParser; % create parser object
addRequired(p,'A');
addRequired(p,'timeA');
addRequired(p,'timeRef');
addParameter(p,'fillin',default_fillin);

parse(p,A,timeA,timeRef,vargs{:})
end

end