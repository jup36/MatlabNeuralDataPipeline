% test the compression scheme

figure(1); colormap(gray);

for i=1:2000
    subplot(211)
    imagesc(seq_image{i}(roi(1,1):roi(1,2),roi(2,1):roi(2,2)),[0 255])
%     imagesc(seq_image{i}(:,:),[0 255])
    subplot(212)
    

    if i==1
        currentImage = first2000compressed.modalImage(roi(1,1):roi(1,2),roi(2,1):roi(2,2));
    else
        deltas             = double(first2000compressed.deltas(i).values);
        indices            = double(first2000compressed.deltas(i).indices);
        tmpImage           = first2000compressed.modalImage(:,:);
        tmpImage(indices)  = tmpImage(indices) + deltas;
        currentImage       = tmpImage(roi(1,1):roi(1,2),roi(2,1):roi(2,2));            
%         currentImage       = tmpImage(:,:);            
    end
    imagesc(currentImage,[0 255]);
%     imagesc(seq_image{i}(roi(1,1):roi(1,2),roi(2,1):roi(2,2))-currentImage,[0 10])
    drawnow;
end