function [fTrj3d, sTrj3d] = stereoTriangulation(fXY, sXY, pr)
% performs triangulation using the stereo-triangulation data (pr) 

[fTrj3d, sTrj3d] = stereo_triangulation(fXY,sXY,pr.om,pr.T,pr.fc_left,pr.cc_left,pr.kc_left,...
pr.alpha_c_left,pr.fc_right,pr.cc_right,pr.kc_right,pr.alpha_c_right);

%plot3(fTrj3d(1,:),fTrj3d(2,:),fTrj3d(3,:))
%plot3(sTrj3d(1,:),sTrj3d(2,:),sTrj3d(3,:))

end
