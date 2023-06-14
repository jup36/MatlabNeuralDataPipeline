% =========================================================
% Read nSamp timepoints from the binary file, starting
% at timepoint offset samp0. The returned array has
% dimensions [nChan,nSamp]. Note that nSamp returned
% is the lesser of: {nSamp, timepoints available}.
%
function dataArray = ReadBin(samp0, nSamp, meta, binName, path)

    nChan = str2double(meta.nSavedChans);

    nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan); 
    % nFileSamp equals to the total number of samples, which equals to time (s) x sampling rate
    % int16 16(2) Bits(Bytes)
    
    samp0 = max(samp0, 0); 
    nSamp = min(nSamp, nFileSamp - samp0); % sampling rate 

    sizeA = [nChan, nSamp]; % # of channels, sampling rate

    fid = fopen(fullfile(path, binName), 'rb');     % fopen opens the file, filename, for binary read access, and returns an integer file identifier equal to or greater than 3. 
    fseek(fid, samp0 * 2 * nChan, 'bof');           % move to specified position in file; fseek sets the file position indicator offset bytes from origin in the specified file. 
    dataArray = fread(fid, sizeA, 'int16=>double'); % readout the binfile with the form of the input changed from int16 to double
    fclose(fid); 
end % ReadBin

