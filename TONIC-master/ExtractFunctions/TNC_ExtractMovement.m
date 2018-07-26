function [outputStructure] = TNC_ExtractMovement(imageArray,roiStruc,invert,compressed)
% FUNCTION DETAILS: Simple function to try and extract movement in a compressed or uncompressed image array loaded into memory.
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: people.janelia.org/dudmanj/html/projects.html
% _________________________________________________________________________

tstart = tic;

    % how many rois for each image?
    numROIs = length(roiStruc);
    
    plotFlag=0;
    
    % set up preallocation
    
    outputStructure.compressed = compressed;

    for indROI = 1:numROIs
        
        clear tmpVect
        
        roi(1,1) = roiStruc(indROI).rowLo;
        roi(1,2) = roiStruc(indROI).rowHi;
        roi(2,1) = roiStruc(indROI).colLo;
        roi(2,2) = roiStruc(indROI).colHi;        
    
        if compressed==1
            dimsOfImage = size(imageArray.modalImage(roi(1,1):roi(1,2),roi(2,1):roi(2,2)));
            totalPixels = dimsOfImage(1,1).*dimsOfImage(1,2);
            totalFrames = size(imageArray.deltas,2);
            outputStructure.varR = zeros(1,totalFrames);
            outputStructure.varC = zeros(1,totalFrames);        
        else
            dimsOfImage = size(imageArray{1,1}(roi(1,1):roi(1,2),roi(2,1):roi(2,2)));
            totalPixels = dimsOfImage(1,1).*dimsOfImage(1,2);
            totalFrames = size(imageArray,2);
            outputStructure.varR = zeros(1,totalFrames);
            outputStructure.varC = zeros(1,totalFrames);
        end

        outputStructure.totalFrames = totalFrames;

        for i = 1:totalFrames
        
            if compressed==1
                if i==1
                    currentImage = imageArray.modalImage(roi(1,1):roi(1,2),roi(2,1):roi(2,2));
                else
                    deltas             = double(imageArray.deltas(i).values);
                    indices            = double(imageArray.deltas(i).indices);
                    tmpImage           = imageArray.modalImage(:,:);
                    if size(deltas)>0
                        tmpImage(indices)  = tmpImage(indices) + deltas;
                    end
                    currentImage       = tmpImage(roi(1,1):roi(1,2),roi(2,1):roi(2,2));            
                end
            else
                currentImage = imageArray{1,i}(roi(1,1):roi(1,2),roi(2,1):roi(2,2));            
            end

            if invert==1
                derivative = 255 - currentImage;       
            else
                derivative = currentImage;
            end

            cumSumRows = cumsum(derivative,1);
            tCumSumR = cumSumRows(size(cumSumRows,1),:);
            cumSumCols = cumsum(derivative,2);
            tCumSumC = cumSumCols(:,size(cumSumRows,2))';

            Rx = 1:length(tCumSumR);
            Rxx = 1:0.1:length(tCumSumR);
            Cx = 1:length(tCumSumC);
            Cxx = 1:0.1:length(tCumSumC);
            tmpSplineR = spline(Rx,tCumSumR,Rxx);
            tmpSplineC = spline(Cx,tCumSumC,Cxx);

            lotCumSumR = decimate(tmpSplineR,30);
            lotCumSumC = decimate(tmpSplineC,30);

            %%%%%%%% SHOULD NOT SAVE THESE TO CONSERVE SPACE        
            %         outputStructure.cumSumRows{1,i}(:,:) = cumSumRows;
                    tmpVect.tCumSumR(i,:) = tCumSumR;
            %         outputStructure.cumSumCols{1,i}(:,:) = cumSumCols;
                    tmpVect.tCumSumC(i,:) = tCumSumC;
            %%%%%%%% SHOULD NOT SAVE THESE TO CONSERVE SPACE        

            outputStructure.roi(indROI).lotCumSumR(i,:) = lotCumSumR;
            outputStructure.roi(indROI).lotCumSumC(i,:) = lotCumSumC;

            if i>1

                outputStructure.roi(indROI).xcorrR(i,:) = xcorr(lotCumSumR,outputStructure.roi(indROI).lotCumSumR(i-1,:),30,'unbiased');
                outputStructure.roi(indROI).xcorrC(i,:) = xcorr(lotCumSumC,outputStructure.roi(indROI).lotCumSumC(i-1,:),30,'unbiased');

                outputStructure.roi(indROI).varR(1,i) = var((tCumSumR-tmpVect.tCumSumR(i-1,:)),0,2);
                outputStructure.roi(indROI).varC(1,i) = var((tCumSumC-tmpVect.tCumSumC(i-1,:)),0,2);
                
                outputStructure.roi(indROI).maxVR(1,i) = max(outputStructure.roi(indROI).xcorrR(i,:));
                outputStructure.roi(indROI).maxLR(1,i) = find(outputStructure.roi(indROI).xcorrR(i,:)==outputStructure.roi(indROI).maxVR(1,i),1);
                outputStructure.roi(indROI).maxVC(1,i) = max(outputStructure.roi(indROI).xcorrC(i,:));
                outputStructure.roi(indROI).maxLC(1,i) = find(outputStructure.roi(indROI).xcorrC(i,:)==outputStructure.roi(indROI).maxVC(1,i),1);

                if plotFlag==1

%                     subplot(3,4,[4,8]);
%                     plot(lotCumSumC-outputStructure.lotCumSumC(i-1,:),1:length(lotCumSumC),'k');
%                     axis([-5000 5000 1 length(lotCumSumC)]);
% 
%                     subplot(3,4,9:11);
%                     plot(1:length(lotCumSumR),lotCumSumR-outputStructure.roi(indROI).lotCumSumR(i-1,:),'r');
%                     axis([1 length(lotCumSumR) -5000 5000]);
% 
%                     subplot(3,4,[1:3,5:7])
%                     plot(1:length(outputStructure.roi(indROI).xcorrC(i,:)),outputStructure.roi(indROI).xcorrC(i,:),'k',1:length(outputStructure.roi(indROI).xcorrR(i,:)),outputStructure.xcorrR(i,:),'r',outputStructure.maxLC(1,i),outputStructure.maxVC(1,i),'ko',outputStructure.maxLR(1,i),outputStructure.maxVR(1,i),'ro');
%                     axis([1 120 -1e7 1e7]);
%                     % imshow(imageArray{1,i}(roi(1,1):roi(1,2),roi(2,1):roi(2,2)),[0 255]); hold on;
% 
%                     subplot(3,4,12);
%                     plot(1:totalFrames,outputStructure.roi(indROI).varR(1,:),'r',1:totalFrames,outputStructure.roi(indROI).varC(1,:),'k');
%                     axis([1 totalFrames 0 1e6]);  
% 
%                     drawnow; hold off;    
%                     pause(0.01);

                else

                    waitbar(i./totalFrames,waith);

                end

            else

                if plotFlag==1
                    figure(1);
                    clf;
                else
                    waith = waitbar(0,'Fraction of frames analyzed...');
                    drawnow;                
                end

            end

        end
        
        if plotFlag==0
            close(waith);
        end
        
    end

    time = toc(tstart);
    disp( sprintf('Data extraction took: %g [sec]',time) );
    