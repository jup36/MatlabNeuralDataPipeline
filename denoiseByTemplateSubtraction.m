function [ denoisedX ] = denoiseByTemplateSubtraction( X, refPoint, varargin  )
%This function receives an (analog) signal X contaminated with an artifact of which 
% occurrences are referenced by refPoints (i.e., one knows exactly when artifacts occurr, 
% e.g., artifact is due to a solenoid activation). The contamination can then be
% eliminated by subtracting the centered artifact from the raw signal at the
% artifact occurrences. 
% If not specified, default pre- and post-paddings are used. 
% One can specify pre- and/or post-paddings. 
% X must be 1-D. 

% check the dimension of X is 1-d 
if ~(size(X,1)==1 || size(X,2)==1)
   error('The input signal has to be 1-d!') 
end

% get prePad and postPad
if nargin == 2
    prePad  = 20; % default pre-padding (# of points before each reference point)
    postPad = 60; % default post-padding (# of points after each reference point)
elseif nargin == 3
    prePad = varargin{1}; % user-input prePad
    postPad = 60;
elseif nargin == 4
    prePad  = varargin{1}; % user-input prePad
    postPad = varargin{2}; % user-input postPad
end

% conduct template subtraction 
denoisedX = X; 
for i = 1:length(refPoint) % increment all the reference points
    if refPoint(i)-prePad>0 && refPoint(i)+postPad<length(X) % if not out-of-bounds
        centeredArtifact = X(refPoint(i)-prePad+1:refPoint(i)+postPad)-mean(X(refPoint(i)-prePad:refPoint(i)-1)); % get centered artifact
        denoisedX(refPoint(i)-prePad+1:refPoint(i)+postPad) = X(refPoint(i)-prePad+1:refPoint(i)+postPad)-centeredArtifact; % subtract artifact from the raw signal to get the denoised signal
    end
end
clearvars i

end
