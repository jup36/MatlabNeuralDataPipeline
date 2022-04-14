function thetaInRadian = angleTwoVecRad(u,v)
    cosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    thetaInDegrees = real(acosd(cosTheta));
    thetaInRadian = deg2rad(thetaInDegrees); 
end