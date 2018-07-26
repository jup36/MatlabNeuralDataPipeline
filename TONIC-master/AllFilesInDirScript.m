% Template to apply a matlab pipeline to a directory of files:
fileExtension = '*.tif'

% Prepare the list all files
totalFilesList      = dir(fileExtension); % creates a structure, use the "name" field
totalFilesCount     = size(totalFilesList,1);

% Create the data structure to hold all of the files

% First print the help of all the high level functions
for j=1:totalFilesCount

    disp(['Analyzing file ' num2str(j) ' of ' num2str(totalFilesCount) ' | ' totalFilesList(j).name]);
    
    currAnalysisWrapper.file(j).name = totalFilesList(j).name;
    currAnalysisWrapper.file(j).date = totalFilesList(j).date;

    currImg = imread(totalFilesList(j).name);

    figure(1);
    imagesc(currImg);
    
    if j==1
        refImg = currImg;
    end
    
    corrScore(j) = corr2(refImg,currImg);
    
    cc = normxcorr2(refImg,currImg); 
    [max_cc, imax] = max(abs(cc(:)));
    
    if j==1
        autocorrPeakLoc = imax;
        [x0,y0] = ind2sub(size(currImg),autocorrPeakLoc)
    end
    
    [x,y] = ind2sub(size(currImg),imax);
    corrShift(j) = pdist([x,y;x0,y0]);
    
end

currAnalysisWrapper.analysis.corr = corrScore;
currAnalysisWrapper.analysis.shift = corrShift;
