function imagescJP( dataMat, colorScheme, colorAxis )
%This is just a helper function for imagesc to simplify the main code 

imagesc(dataMat);

set(gca,'TickDir','out')
axis tight
caxis(colorAxis);
axis xy;
colormap(colorScheme);

end

