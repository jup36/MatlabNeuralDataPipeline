arrayType       = 'NN_b64';
dispOn          = 0;
chunk           = 180;
basePath        = '/Volumes/dudmanlab/users/Babita/Head Fixed/';
[num,txt,raw]   = xlsread('~/Desktop/Remaining-Striatal-Files.xls');
numFiles        = size(txt,1)

for q = 74:numFiles

    baseName    = cell2mat(txt(q));
    path        = [basePath baseName(1:3) '/Neural Recordings/'];
    fileNames   = dir([path baseName '*.ns5']);
    targetName  = fileNames(1).name;
    fileNameStr = [path targetName];
    
    [ContData] = TNC_MoverBehaviorExtract(fileNameStr,['/Volumes/dudmanlab/users/Babita/ToAnalyze/' baseName]);
    clear ContData

    [sessionStruct,featStruct] = TNC_SSP_Features(fileNameStr,arrayType,chunk,dispOn,['/Volumes/dudmanlab/users/Babita/ToAnalyze/' baseName]);    
    clear sessionStruct featStruct;
    
end 