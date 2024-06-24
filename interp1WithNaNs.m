function newDatInterp = interp1WithNaNs(origTs, origDat, newTs)
% origTs = faLickTs(faLickFrameI); 
% origDat = raw_pupil_area_lickAlign; 
% newTs = lickPeth; 

if sum(isnan(origDat))>0
    if sum(isnan(origDat))/length(origTs)<0.2
        nanI = isnan(origDat); 
        datWoNaN = interp1(origTs(~nanI), origDat(~nanI), origTs, 'linear', 'extrap'); % fill in NaNs
    else
        datWoNaN = origDat; % do not fill in for NaNs 
    end
else
    datWoNaN = origDat; 
end

if ~isequal(origTs, newTs)
   newDatInterp = interp1(origTs, datWoNaN, newTs, 'linear', 'extrap'); % fill in NaNs
else
    newDatInterp = datWoNaN; 
end


end