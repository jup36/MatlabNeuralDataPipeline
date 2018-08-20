%This script exemplifies how one can get all columnwise combinations of two
% matrices and get the product of all column combinations. 

a = [1     0     0
     1     0     0
     1     0     0
     1     0     0
     0     1     0
     0     1     0
     0     1     0
     0     1     0
     0     0     1
     0     0     1
     0     0     1
     0     0     1];
 b= [0     0     1
     0     1     0
     1     0     0
     0     0     1
     0     1     0
     1     0     0
     0     0     1
     0     1     0
     1     0     0
     0     0     1
     0     1     0
     1     0     0];
 
 out=reshape(repmat(a,size(b,2),1),size(a,1),[]) & repmat(b,1,size(a,2)); 

 