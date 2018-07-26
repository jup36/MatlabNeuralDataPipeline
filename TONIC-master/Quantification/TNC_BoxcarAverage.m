function [boxCarAvg] = TNC_BoxcarAverage(window,data)

halfwin = (window-1)./2;
if size(data,1)>size(data,2)
    tmpData = data';
else
    tmpData = data;
end

boxCarAvg = zeros(1,size(tmpData,2));

for i=1:size(tmpData,2)
    
    if i < halfwin+1
        boxCarAvg(1,i) = mean(tmpData(1,1:i+halfwin),2);
    elseif i+halfwin > size(tmpData,2)
        boxCarAvg(1,i) = mean(tmpData(1,i-halfwin:size(tmpData,2)),2);
    else
        boxCarAvg(1,i) = mean(tmpData(1,i-halfwin:i+halfwin),2);
    end
    
end
    
    