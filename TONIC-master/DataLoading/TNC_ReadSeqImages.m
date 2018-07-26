function seq_image = TNC_ReadSeqImages(seq_info, fid, frames)
% seq_im = read_seq_images(seq_info, fid, frames)
% reads the frames specified by frames, either a range, or -1 for all of
% them, and returns in the seq_im structure

% tstart = tic;

disp(sprintf('Loading %g frames...',length(frames)));

if frames == -1  % if -1, then set frames to entire sequence
    frames = 1:seq_info.NumberFrames;
end

if length(frames) == 1 % if just one frame, do not return a structure
    image_address = 1024 + (frames-1)*seq_info.TrueImageSize;
    status = fseek(fid, image_address, 'bof');
    seq_image = fread(fid, [seq_info.Width, seq_info.Height], 'uint8')';
else
    for j = 1:length(frames)
        image_address = 1024 + (frames(1,j)-1).*seq_info.TrueImageSize;
        status = fseek(fid, image_address, 'bof');
%         seq_image.frame(j).data = fread(fid, [seq_info.Width, seq_info.Height], 'uint8')';
        seq_image{1,j}(:,:) = fread(fid, [seq_info.Width, seq_info.Height], 'uint8')';
    end    
end

% toc(tstart)

