
function p = targetintoreftimePrep( A,timeA,timeRef, vargs )
% parse input, and extract name-value pairs
default_fillin = 0; % by default do not re-read the raw bin file, if done already

p = inputParser; % create parser object
addRequired(p,'A');
addRequired(p,'timeA');
addRequired(p,'timeRef');
addParameter(p,'fillin',default_fillin);

parse(p,A,timeA,timeRef,vargs{:})
end