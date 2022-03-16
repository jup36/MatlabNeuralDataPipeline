function reachAngle = computeReachAngleC( hTrjC, jsXYposC )
%Angle is calculated relative to the straight line projected from the
% initial hand position.  
    distToJs = cellfun(@(a, b) sum((a-repmat(b,1,size(a,2))).^2,1), hTrjC, jsXYposC, 'un', 0);  % compute the distance from the jsXY
    for j = 1:length(distToJs)
        [~,tmpMinI] = min(distToJs{j});
        u=hTrjC{j}(:,tmpMinI)-hTrjC{j}(:,1); % reach vector
        %v=[hTrjC{j}(1,tmpMinI);hTrjC{j}(2,1)]-hTrjC{j}(:,1); % reference vector
        v=[hTrjC{j}(1,1);hTrjC{j}(2,tmpMinI)]-hTrjC{j}(:,1); % reference vector
        reachAngle(j) = angleTwoVectors(u,v); % get the angle between reach and reference vectors 
    end
end
    