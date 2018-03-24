function J = image2sum(I)
% function J = image2sum(I)
%
% I is a 2D-matrix (an image)
%
% useful preprocessing for image normalization
%
% J is the same size as I
%
% J(j,k) = sum(sum(I(1:j,1:k)));
%
% notice that 
% sum(sum(I(j:m,k:n))) = J(m,n)-J(j-1,n)-J(m,k-1)+J(j-1,k-1)
%
% see also image2sumF.m

% cast to double incase I is an integer image
J = double(I);

[r,c] = size(J);

% process the first column
for m = 2:r
    J(m,1) = J(m,1) + J(m-1,1);
end

% process subsequent columns
for n = 2:c
    Isum = 0;
    for m = 1:r
        % keep track of the sum running down this column
        Isum = Isum + J(m,n);
        % add this running sum to everything just to the left
        J(m,n) = J(m,n-1)+Isum;
    end
end
