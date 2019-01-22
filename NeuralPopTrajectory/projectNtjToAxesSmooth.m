function nTjOut = projectNtjToAxesSmooth( nTj, projVecs, gaussianKernel )
%This function takes a neural trajectories in a cell array, and project
% them on the axes defined by the projVecs, and then smooth the trajectories
% using the gaussian kernel provided.

if iscell(nTj)
    if cellfun(@(x) size(x,1), nTj) == size(projVecs,1) % check if the dimensions of nTj and projVecs match!
        for f = 1:length(nTj)
            tmpNtj = projVecs'*nanmean(nTj{f},3); % project the trial-averaged population activity of each fold to the projection matrix
            nTjOut{1,f} = cell2mat(arrayfun(@(ROWIDX) conv(tmpNtj(ROWIDX,:),gaussianKernel,'same'), (1:size(tmpNtj,1))', 'UniformOutput', false)); % smooth
        end
        clearvars f
    else
        error('The dimensions between nTj and projVecs do NOT match!')
    end
else
    error('nTj must be provided as a cell array!')
end

end

