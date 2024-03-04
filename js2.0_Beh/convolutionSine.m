function [timeseries_conv, timeseries_conv_C_meanSub] = convolutionSine(timeseries_org, freqList, sampleRate)
%This function conducts convolution between the original timeseries ('timeseries_org')
% and sine waves of frequencies specified by the 'freqList'.
dur = 0.01; % duration of the sine waves (10 ms). Note that one must use a small window for temporal precision, the expected temporal precision with 10-ms window is about 5 ms.
assert(unique(dur./(1./freqList)>1)); % ensure that the dur is longer than each frequency component's full cycle, in other words, for convolution high-enough frequencies must be used for the 10-ms window
t = 0:1/sampleRate:dur-1/sampleRate;

timeseries_conv_C = cell(length(freqList), 1);
for fq = 1:length(freqList)
    sineW = sin(2*pi*freqList(fq)*t); % build sine wave for convolution
    timeseries_conv_C{fq, 1} = conv(timeseries_org, sineW, 'same');
end

timeseries_conv_C_meanSub = cellfun(@(a) a-nanmean(a), timeseries_conv_C, 'UniformOutput', false);

timeseries_conv = sum(cell2mat(timeseries_conv_C_meanSub), 1);

end