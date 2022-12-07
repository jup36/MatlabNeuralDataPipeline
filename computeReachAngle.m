function reachAngle = computeReachAngle( hTrjC, jsXYpos )
%Angle is calculated relative to the straight line projected from the
% initial hand position.  
    distToJs = cellfun(@(a) sum((a-repmat(jsXYpos,1,size(a,2))).^2,1), hTrjC, 'un', 0);  % compute the distance from the jsXY
    sign = [1,-1]; % take rightward move as negative angles
    for j = 1:length(distToJs)
        [~,tmpMinI] = min(distToJs{j});
        u=hTrjC{j}(:,tmpMinI)-hTrjC{j}(:,1); % reach vector
        %v=[hTrjC{j}(1,tmpMinI);hTrjC{j}(2,1)]-hTrjC{j}(:,1); % reference vector
        v=[hTrjC{j}(1,1);hTrjC{j}(2,tmpMinI)]-hTrjC{j}(:,1); % reference vector
        reachAngle(j) = sign((u(1)>0)+1) * angleTwoVectors(u,v); % get the angle between reach and reference vectors      
    end
    
function thetaInDegrees = angleTwoVectors(u,v)
    cosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    thetaInDegrees = real(acosd(cosTheta));
end    
end
    