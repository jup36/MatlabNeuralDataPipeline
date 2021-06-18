function [hTrj3dCout] = hTrjMedianSubtractRotate(hTrj3dC)
%This function takes the raw hand trajectories and conduct median
% subtraction and rotation
%hTrj3dC = {jkvt(:).hTrjF}; % hand trajectory (original)
if ~exist('cwRot3d','var') % get the calibrated rotation matrix if needed
    load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/collectData','clockwiseRotationMatrix'),'cwRot3d') % load the rotation matrix to re-align hand coordinates
end
fillTr = ~cell2mat(cellfun(@isempty,hTrj3dC,'un',0)); % trials with hand trajectories
medHInit = nanmedian(cell2mat(cellfun(@(a) a(:,1), hTrj3dC(fillTr), 'un', 0)),2); % median initial hand position for subtraction
hTrjMS = cellfun(@(a) a-repmat(medHInit,[1,size(a,2)]), hTrj3dC(fillTr), 'un', 0); % subtract the median initial hand position
hTrj3dC(fillTr) = deal(hTrjMS); % put the median subtracted trajectories in
hTrjR = deal(cellfun(@(a) cwRot3d*a, hTrj3dC(fillTr), 'un', 0)); % clockwise rotation
hTrj3dC(fillTr) = deal(hTrjR); % put rotated trajectories in
hTrj3dCout = hTrj3dC;
end
