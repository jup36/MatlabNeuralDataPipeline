function [alignmentMetrics] = TNC_ExtractImageAlignment(imageArray,invert,compressed)

wholeImage=0;

tstart = tic;

    outputStructure.compressed = compressed;

    if compressed==1
        dimsOfImage = size(imageArray.modalImage);
        totalPixels = dimsOfImage(1,1).*dimsOfImage(1,2);
        totalFrames = size(imageArray.deltas,2);
    else
        dimsOfImage = size(imageArray{1,1}(:,:));
        totalPixels = dimsOfImage(1,1).*dimsOfImage(1,2);
        totalFrames = size(imageArray,2);
    end

    profile50 = zeros(totalFrames,dimsOfImage(1,2));
    corrScore = zeros(1,totalFrames);
    corrShift = zeros(1,totalFrames);
    corrShiftRatio = zeros(1,totalFrames);
    
    outputStructure.totalFrames = totalFrames;
    
    profile50loc = round(dimsOfImage(1,1).*0.5)

    for i = 1:totalFrames
        
        if compressed==1
            if i==1
                currentImage = imageArray.modalImage(:,:);
            else
                deltas             = imageArray.deltas{1,i}(:,:);
                indices            = imageArray.indices{1,i}(:,:);
                tmpImage           = imageArray.modalImage(:,:);
                tmpImage(indices)  = tmpImage(indices) + deltas;
                currentImage       = tmpImage(:,:);
            end
        else
            currentImage = imageArray{1,i}(:,:);            
        end

        if invert==1
            derivative = 255 - currentImage;       
        else
            derivative = currentImage;
        end

        if i==1
            refImg = derivative;
        end        
        
        profile50(i,:) = derivative(profile50loc,:);

        corrScore(i) = corr2(refImg,derivative);

        if wholeImage==1
%             cc = normxcorr2(refImg,derivative); % not meaningful unless different sizes should use xcorr2
            cc = xcorr2(refImg,derivative);
            [max_cc, imax] = max(abs(cc(:)));
            if i==1
                autocorrPeakLoc = imax;
                [x0,y0] = ind2sub(size(derivative),autocorrPeakLoc)
            end

            [x,y] = ind2sub(size(derivative),imax);
            corrShift(i) = pdist([x,y;x0,y0]);
        
        else
            tmp = xcorr(profile50(1,:),profile50(i,:));
            maxCorr = find(tmp==max(tmp));
            corrShift(i) = maxCorr - dimsOfImage(1,2);
            corrShiftRatio(i) = (tmp(maxCorr) - ((tmp(maxCorr-1) + tmp(maxCorr+1))./2))./tmp(maxCorr);
        end
        


%         disp(i);
%         figure(1);
%         subplot(212);
%         plot(1:i,corrShift(1:i),'k',1:i,corrScore(1:i),'r');
%         axis([0 totalFrames -1 4]);
%         subplot(211);
%         plot(1:dimsOfImage(1,2),profile25(i,:),'r',1:dimsOfImage(1,2),profile75(i,:),'k');
%         axis([0 dimsOfImage(1,2) 0 255]);

    end
    
    alignmentMetrics.corrShift = corrShift;
    alignmentMetrics.corrScore = corrScore; 
    alignmentMetrics.corrShiftRatio = corrShiftRatio; 
    
    alignmentMetrics.profiles.profile50 = profile50;
    
    time = toc(tstart);

disp( sprintf('Data extraction took: %g [sec]',time) );

% recap
% 
% for i=1:1000
%     subplot(211);
%     imagesc(imageArray{1,i}(:,:));
%     subplot(212);
%     plot(alignmentMetrics.profiles.profile50(i,:));
%     drawnow;
% end

figure(2)
subplot(311);
plot(alignmentMetrics.corrScore);
subplot(312)
plot(alignmentMetrics.corrShift);
subplot(313)
plot(alignmentMetrics.corrShiftRatio);


