function subdirinfo = dirsub(filePath, fileName)

dirinfo = dir(filePath);
dirinfo(~[dirinfo.isdir]) = [];

subdirinfo = {};
for K = 1 : length(dirinfo)
    thisdir = fullfile(dirinfo(K).folder,dirinfo(K).name);
    subdirinfo = [subdirinfo; dir(fullfile(thisdir, fileName))]; 
end

end