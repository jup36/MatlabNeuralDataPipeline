function thresh = calcThresh( currSamples, qqFactor, nBit )

    quirogaDenom = 0.6745;
    
    currDev = abs(currSamples);
    
    nEdge = nBit;
    hist_edges = 0:nEdge;
    [counts,~] = histcounts(currDev,hist_edges);
    nValues = sum(counts);
    
    medFound = 0;
    iE = 1;
    while( ~medFound && iE < nEdge )
        sumLow = sum(counts(1:iE));
        sumHigh = sum(counts(iE+1:nEdge));
        if( sumHigh < sumLow )
            medFound = 1;
            currMed = hist_edges(iE);
            %calculate the estimated median within the median bin using a
            %linear interpolation between the lower and upper edge of the
            %median bin. This is known as estimating the mean for grouped
            %data -- the groups in this case are the low bits of the int16
            %values.
            
            lb = currMed; %lower bound of bin containing median
            if iE > 1
                yDist = nValues/2 - sum(counts(1:iE-1));
            else
                yDist = nValues/2; %if median is in zero bin, sum of bins below med = 0;
            end
            estMed = lb + yDist/counts(iE);
            
            %fprintf( '%d\t%d\t%.2f\n', selectChan(i), currMed, estMed(selectChan(i)) );
        else
            iE = iE + 1;
        end
    end
    
    if( medFound == 0 ) 
            fprintf( 'median not found. problem reading data?\n' );
            thresh = -1;
            return;
    else
        thresh = qqFactor*(estMed/quirogaDenom);
    end

end