function [] = TNC_ExportRasterToIgor(rasterStruct,tWindow,path) 

pathname = [path '.h5']
size(dataMat,2)

% numTrials = size(PopData.session(sessInd).unit(unitInd).rLITE.raster.trial,2);
numTrials = size(PopData.session(sessInd).unit(unitInd).rTONE.raster.trial,2);

for n = 1:numTrials
%     thisStamps = PopData.session(sessInd).unit(unitInd).rLITE.raster.trial(n).ts;
    thisStamps = PopData.session(sessInd).unit(unitInd).rTONE.raster.trial(n).ts;

    spkSTAMPS = [spkSTAMPS thisStamps'];
    trialLIST = [trialLIST ones(1,numel(thisStamps)).*n];
    
end

exportToIgor = [trialLIST' spkSTAMPS'];
