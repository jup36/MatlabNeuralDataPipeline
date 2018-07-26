function [] = TNC_PlotContour(c,lwidth,currColor,gradation,figNum,clear,zOffset)


figure(figNum); 
if clear
    clf; 
end
hold on;

numElem = size(c,2);

i = 1;
j = 2;
count = 1;
currPathLength = c(2,i);
currLevel = c(1,i);

while 1

    currPathLength = c(2,i);
    currLevel = c(1,i);

    if gradation~=0
        currL = lwidth.*(currLevel./gradation);
    else
        currL = lwidth;
    end

    if zOffset>0
        plot3(c(1,j:j+currPathLength-1),c(2,j:j+currPathLength-1),ones(1,currPathLength).*zOffset,'Color',currColor,'LineWidth',currL);        
    else
        plot(c(1,j:j+currPathLength-1),c(2,j:j+currPathLength-1),'Color',currColor,'LineWidth',currL);
    end

    
    count = count + 1;

    i = j+currPathLength;
    j = i+1;
    
    
    if i>numElem
        %disp('complete.');
        break
    end
end