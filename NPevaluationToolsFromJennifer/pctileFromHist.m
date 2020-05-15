function [histPctile] = pctileFromHist( hData, p )

%calculate the bin corresponding to a given percentile p (value 0 to 1)
%from a histogram

nBin = numel(hData);

if sum(hData) > 0
    targVal = p*sum(hData);

    i = 1;
    currVal = 0;

    while( currVal <= targVal && i < nBin )
        currVal = currVal + hData(i);
        i = i + 1;
    end

    histPctile = int32(i);
else
    histPctile = 1;
end
  
end