function imagescWithSmoothing( dataMat, colorScheme, colorAxis, smoothFactor )
%This is just a helper function for imagesc to simplify the main code 

% smoothing
if smoothFactor > 1
   smDataMat = zeros(size(dataMat)); 
   for i = 1:size(dataMat,1)
        smDataMat(i,:) = smooth(dataMat(i,:),smoothFactor); 
   end
else
    smDataMat = dataMat;
end

imagesc(smDataMat);

set(gca,'TickDir','out')
axis tight
caxis(colorAxis);
axis xy;
colormap(colorScheme);

end

