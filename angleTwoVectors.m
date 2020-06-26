function thetaInDegrees = angleTwoVectors(u,v)
    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    thetaInDegrees = real(acosd(CosTheta));
end