function imagescJP( dataMat, colorScheme, colorAxis, varargin )
%This is just a helper function for imagesc to simplify the main code 

p = parse_input_imagescJP( dataMat, colorScheme, colorAxis, varargin );

imagesc( p.Results.X, p.Results.Y, dataMat);

set(gca,'TickDir','out')
axis tight
caxis(colorAxis);
axis xy;
colormap(colorScheme);
if p.Results.cBar
    colorbar
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED HELPER FUNCTIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = parse_input_imagescJP( filePath, fileName, varName, vargs ) 

default_X = [1,size(dataMat,2)]; % define X corners for imagesc
default_Y = [1,size(dataMat,1)]; % define Y corners for imagesc
default_cBar = true;

p = inputParser; % create parser object
addRequired(p,'dataMat');
addRequired(p,'colorScheme');
addRequired(p,'colorAxis');

addParameter(p,'X', default_X)
addParameter(p,'Y', default_Y)
addParameter(p,'cBar', default_cBar)

parse(p, filePath, fileName, varName, vargs{:})

end

end

