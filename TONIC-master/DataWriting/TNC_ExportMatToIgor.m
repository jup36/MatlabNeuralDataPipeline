function [] = TNC_ExportMatToIgor(dataMat,tWindow,path) 

pathname = [path '.h5']
size(dataMat,2)

for i = 1:size(dataMat,1)
    
    currData    = dataMat(i,tWindow);

    if i==1
        name = sprintf('/RecordA%g',i-1);
        hdf5write(pathname,name,currData);
    else
        name = sprintf('/RecordA%g',i-1);
        hdf5write(pathname,name,currData,'WriteMode','append');
    end 
        
end