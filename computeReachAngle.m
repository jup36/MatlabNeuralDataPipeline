function reachAngle = computeReachAngle( hTrjC, jsXYpos )
    distToJs = cellfun(@(a) sum((a-repmat(jsXYpos,1,size(a,2))).^2,1), hTrjC, 'un', 0);  % compute the distance from the jsXY
    for j = 1:length(distToJs)
        [~,tmpMinI] = min(distToJs{j});
        u=hTrjC{j}(:,tmpMinI)-hTrjC{j}(:,1); % reach vector
        v=[hTrjC{j}(1,tmpMinI);hTrjC{j}(2,1)]-hTrjC{j}(:,1); % reference vector
        reachAngle(j) = angleTwoVectors(u,v); % get the angle between reach and reference vectors 
    end
end
    