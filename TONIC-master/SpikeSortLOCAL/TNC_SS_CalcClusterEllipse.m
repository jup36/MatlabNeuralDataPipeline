function [eX,eY] = TNC_SS_CalcClusterEllipse(cntX,stdX,cntY,stdY,numStd,angle,scaler)

    steps=36;

    beta = -angle * (pi / 180);
    sinbeta = sin(beta);
    cosbeta = cos(beta);

    alpha = linspace(0, 360, steps)' .* (pi / 180);
    sinalpha = sin(alpha);
    cosalpha = cos(alpha);

    eX = cntX + ((stdX.*numStd) * cosalpha * cosbeta - (stdY.*numStd) * sinalpha * sinbeta);
    eY = cntY + ((stdX.*numStd) * cosalpha * sinbeta + (stdY.*numStd) * sinalpha * cosbeta);
    
    %% Include a rule akin to lognormal like distributions
    
    distFrom0 = sqrt(eX.^2+eY.^2);
    normDist = (distFrom0-min(distFrom0)) ./ (max(distFrom0)-min(distFrom0));
    
    posX =   find(eX>0);
    posY =   find(eY>0);
    
    for j=1:numel(eX)
        if eX(j)-cntX>0
            eX(j) = eX(j) + ((eX(j)-cntX).*scaler.*normDist(j));
        else
            eX(j) = eX(j) - (abs(eX(j)-cntX).*scaler.*normDist(j));            
        end
        
        if eY(j)-cntY>0
            eY(j) = eY(j) + ((eY(j)-cntY).*scaler.*normDist(j));
        else
            eY(j) = eY(j) - (abs(eY(j)-cntY).*scaler.*normDist(j));            
        end
        
    end
    
