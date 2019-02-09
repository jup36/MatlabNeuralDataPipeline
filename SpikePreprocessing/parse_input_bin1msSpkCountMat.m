function p = parse_input_bin1msSpkCountMat( SpkCountMat1ms, binSize, stepSize, vargs )
%parse input, and extract name-value pairs for the main function 'bin1msSpkCountMat.m'

default_align = 'center';

p = inputParser; % create parser object

addRequired(p,'SpkCountMat1ms'); % SpkCountMat to be binned
addRequired(p,'binSize');  % binSize
addRequired(p,'stepSize'); % stepSize

addParameter(p,'align', default_align) % either 'center' to take preceding and following bins or 'subsequent' to just take following bins

parse(p, SpkCountMat1ms, binSize, stepSize, vargs{:})
end