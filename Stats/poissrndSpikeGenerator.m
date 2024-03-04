function yhatUnitByTimeRepPoiss = poissrndSpikeGenerator(yhatUnitByTime)
    % replace negative values in yhatUnitByTime to zeros.  
    yhatUnitByTime(yhatUnitByTime<0)=0; 
    % poiss random sampling
    yhatUnitByTimeRepPoiss = poissrnd(yhatUnitByTime); 
end