% write out some hdf5 data for the whole data array
% 
% reMap = ...
% [5,13,21,29,37,45,53,61;...
%  4,12,20,28,36,44,52,60;...
%  6,14,22,30,38,46,54,62;...
%  3,11,19,27,35,43,51,59;...
%  7,15,23,31,39,47,55,63;...
%  2,10,18,26,34,42,50,58;...
%  8,16,24,32,40,48,56,64;...
%  1,9,17,25,33,41,49,57];

% tmp = data.Data(reMap,:);
% tmp = data.Data;
% tmp = dataAinp.Data;
clear currDat*

for i = 1:663
    
%     currData = allCSaligned.sorted(:,i);
    currData    = allCSaligned.sort.susPeak.sortedZ(:,i);
    currData2   = allCSEXaligned.sort.susPeak.sortedZ(:,i);

    if i==1
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../testExportAN1.h5',name,currData);

        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../testExportAN2.h5',name,currData2);
    else
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../testExportAN1.h5',name,currData,'WriteMode','append');
        name = sprintf('/RecordA%g',i-1);
        hdf5write('../../testExportAN2.h5',name,currData2,'WriteMode','append');
    end 
        
end

    % also write out the cluster indices:
    name = sprintf('/Ids%g',i-1);
    hdf5write('../../testExportAN1.h5',name,allIds,'WriteMode','append');
