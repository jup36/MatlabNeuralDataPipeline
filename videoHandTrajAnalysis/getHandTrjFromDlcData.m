function getHandTrjFromDlcData(filePath, hTrjPath)

%filePath = '/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081419';
%filePath = 'S:\Junchol_Data\JS2p0\WR40_081419'; % data folder (beefcake)
%hTrjPath = '/Volumes/dudmanlab/dilucid-drop-output/mouseJoystick3/081419_WR40';
%hTrjPath = 'T:\dilucid-drop-output\mouseJoystick3\081419_WR40'; % hand trajectory folder

%% load streo-calibration data 
stereoCalibFile = 'C:\Users\parkj\Documents\TOOLBOX_calib\Calib_Results_stereo_011619.mat'; % caltech camera calibration toolbox
if exist(fullfile(stereoCalibFile),'file')==2
    calibParams = load(fullfile(stereoCalibFile),'om','T','fc_left','cc_left','kc_left','alpha_c_left','fc_right','cc_right','kc_right','alpha_c_right');
else
    disp('Point to "Calib_Results_stereo_011619.mat" with stereo camera calibration data!')
    [calibFileSelect,calibPathSelect] = uigetfile('C:\Users\parkj\Documents\TOOLBOX_calib');
    calibParams = load(fullfile(calibPathSelect,calibFileSelect),'om','T','fc_left','cc_left','kc_left','alpha_c_left','fc_right','cc_right','kc_right','alpha_c_right'); % load stereo camera calibration data
end
    
%% collect jsTime1k_Kinematics data
jsKinFile = dir(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles.mat')); 

if length(jsKinFile)==1
    load(fullfile(filePath,'jsTime1k_Kinematics_VideoFiles.mat'),'jsTime1k_KV'); % load jsTime1k_KV
else
    disp('Point to "jsTime1k_Kinematics_VideoFiles.mat" with the variable "jsTime1k_KV"!')
    [jsKinFileSelect,jsKinPathSelect] = uigetfile(filePath);
    load(fullfile(jsKinPathSelect,jsKinFileSelect),'jsTime1k_KV'); % load jsTime1k_KV
end

fVideoInfo = {jsTime1k_KV(:).fVideoInfo}';
sVideoInfo = {jsTime1k_KV(:).sVideoInfo}';
vFrameTime = {jsTime1k_KV(:).vFrameTime}';
vUseFrameIdx = {jsTime1k_KV(:).vUseFrameIdx}';
clearvars jsTime1k_KV

%% get the list of dlc output csv files
if exist(fullfile(hTrjPath),'dir')==7
    allCsv = dir(fullfile(hTrjPath,'*.csv')); % all trial file
else
    disp('Point to the folder with dlc output csv files (E.g. "cam0_cam_0_2019_08_14_14_19_01.csv")!')
    hTrjPath = uigetdir('T:\dilucid-drop-output\mouseJoystick3');
    allCsv = dir(fullfile(hTrjPath,'*.csv')); % all trial file
end
    
%% get frame dimensions
fronIdx = find(cell2mat(cellfun(@isstruct, fVideoInfo,'Un',0))==1,1,'first'); % get the 1st valid fron video file
sideIdx = find(cell2mat(cellfun(@isstruct, sVideoInfo,'Un',0))==1,1,'first'); % get the 1st valid side video file

[~,vFolderName] = fileparts(hTrjPath);
fVid = VideoReader(fullfile(filePath,vFolderName,fVideoInfo{fronIdx}.name)); % read the 1st front video to get frame info
sVid = VideoReader(fullfile(filePath,vFolderName,sVideoInfo{sideIdx}.name)); % read the 1st side video to get frame info
wd = fVid.Width;  % frame width
ht = fVid.Height; % frame height

%% Read pixel values of each body part from CSV files
fronCnt = 0; % count front cam files
sideCnt = 0; % count side cam files
for f = 1:length(allCsv)
    tmpTrj.FileName = fullfile(allCsv(f).folder,allCsv(f).name);
    tmpCsv = readtable(tmpTrj.FileName);
    tmpCsvCell = table2cell(tmpCsv);
    
    if f == 1
        % locate variables in the csv file
        fg1Idx = cellfun(@(a) strcmpi(a,'Finger1'), tmpCsvCell(1,:));
        fg2Idx = cellfun(@(a) strcmpi(a,'Finger2'), tmpCsvCell(1,:));
        hdIdx = cellfun(@(a) strcmpi(a,'hand'), tmpCsvCell(1,:));
        jsTIdx = cellfun(@(a) strcmpi(a,'Joystick1'), tmpCsvCell(1,:));
        jsBIdx = cellfun(@(a) strcmpi(a,'Joystick2'), tmpCsvCell(1,:));
        
        xIdx = cellfun(@(a) strcmpi(a,'x'), tmpCsvCell(2,:));
        yIdx = cellfun(@(a) strcmpi(a,'y'), tmpCsvCell(2,:));
        lIdx = cellfun(@(a) strcmpi(a,'likelihood'), tmpCsvCell(2,:));
    end
    
    tmpTrj.fg1X = cellfun(@str2double, tmpCsvCell(3:end,fg1Idx & xIdx)); % finger1 x
    tmpTrj.fg1Y = cellfun(@str2double, tmpCsvCell(3:end,fg1Idx & yIdx)); % finger1 y
    tmpTrj.fg1lh = cellfun(@str2double, tmpCsvCell(3:end,fg1Idx & lIdx)); % finger1 likelihood
    tmpTrj.fg1IntX = interplh(tmpTrj.fg1X, tmpTrj.fg1lh, wd, false); % interpolate based on data likelihood (confidence)
    tmpTrj.fg1IntY = interplh(tmpTrj.fg1Y, tmpTrj.fg1lh, ht, false); % interpolate based on data likelihood (confidence)
    % trj = tmpTrj.fg1Y; trjlh = tmpTrj.fg1lh; lm = ht;
    
    tmpTrj.fg2X = cellfun(@str2double, tmpCsvCell(3:end,fg2Idx & xIdx)); % finger2 x
    tmpTrj.fg2Y = cellfun(@str2double, tmpCsvCell(3:end,fg2Idx & yIdx)); % finger2 y
    tmpTrj.fg2lh = cellfun(@str2double, tmpCsvCell(3:end,fg2Idx & lIdx)); % finger2 likelihood
    tmpTrj.fg2IntX = interplh(tmpTrj.fg2X, tmpTrj.fg2lh, wd, false); % interpolate based on data likelihood (confidence)
    tmpTrj.fg2IntY = interplh(tmpTrj.fg2Y, tmpTrj.fg2lh, ht, false); % interpolate based on data likelihood (confidence)
    
    tmpTrj.hdX = cellfun(@str2double, tmpCsvCell(3:end,hdIdx & xIdx)); % hand x
    tmpTrj.hdY = cellfun(@str2double, tmpCsvCell(3:end,hdIdx & yIdx)); % hand y
    tmpTrj.hdlh = cellfun(@str2double, tmpCsvCell(3:end,hdIdx & lIdx)); % hand likelihood
    tmpTrj.hdIntX = interplh(tmpTrj.hdX, tmpTrj.hdlh, wd, false); % interpolate based on data likelihood (confidence)
    tmpTrj.hdIntY = interplh(tmpTrj.hdY, tmpTrj.hdlh, ht, false); % interpolate based on data likelihood (confidence)
    % trj = tmpTrj.hdY; trjlh = tmpTrj.hdlh; lm = ht;
    % trj = tmpTrj.hdX; trjlh = tmpTrj.hdlh; lm = wd;
    
    % jsT: joystick Top
    tmpTrj.jsTX = cellfun(@str2double, tmpCsvCell(3:end,jsTIdx & xIdx)); % finger1 x
    tmpTrj.jsTY = cellfun(@str2double, tmpCsvCell(3:end,jsTIdx & yIdx)); % finger1 y
    tmpTrj.jsTlh = cellfun(@str2double, tmpCsvCell(3:end,jsTIdx & lIdx)); % finger1 likelihood
    [tmpTrj.jsTIntX,tmpTrj.jsTfstPtX] = interplh(tmpTrj.jsTX, tmpTrj.jsTlh, wd, true); % interpolate based on data likelihood (confidence)
    [tmpTrj.jsTIntY,tmpTrj.jsTfstPtY] = interplh(tmpTrj.jsTY, tmpTrj.jsTlh, ht, true); % interpolate based on data likelihood (confidence)
    % trj = tmpTrj.jsTY; trjlh = tmpTrj.jsTlh; lm = ht;
    % trj = tmpTrj.jsTX; trjlh = tmpTrj.jsTlh; lm = wd;
    
    % jsB: joystick Bottom
    tmpTrj.jsBX = cellfun(@str2double, tmpCsvCell(3:end,jsBIdx & xIdx)); % finger1 x
    tmpTrj.jsBY = cellfun(@str2double, tmpCsvCell(3:end,jsBIdx & yIdx)); % finger1 y
    tmpTrj.jsBlh = cellfun(@str2double, tmpCsvCell(3:end,jsBIdx & lIdx)); % finger1 likelihood
    [tmpTrj.jsBIntX,tmpTrj.jsBfstPtX] = interplh(tmpTrj.jsBX, tmpTrj.jsBlh, wd, true); % interpolate based on data likelihood (confidence)
    [tmpTrj.jsBIntY,tmpTrj.jsBfstPtY] = interplh(tmpTrj.jsBY, tmpTrj.jsBlh, ht, true); % interpolate based on data likelihood (confidence)
    
    if contains(tmpTrj.FileName,'cam0','IgnoreCase',true)
        fronCnt = fronCnt + 1;
        fron(fronCnt) = tmpTrj;
    elseif contains(tmpTrj.FileName,'cam1','IgnoreCase',true)
        sideCnt = sideCnt + 1;
        side(sideCnt) = tmpTrj;
    end
    
    fprintf('processed file #%d\n', f);
    %plot(tmpTrj.hdIntX, tmpTrj.hdIntY)
end
save(fullfile(filePath,'rawFHtrj'), 'fron','side'); % save the raw pixel values
% load(fullfile(filePath,'rawFHtrj'), 'fron','side'); % save the raw pixel values

%% get 3-d trajectories with streoTriangulation and assign them to corresponding trials in kv
fList = {fron(:).FileName}; % front video list
sList = {side(:).FileName}; % side video list

for t = 1:length(fVideoInfo) % increment trials of jsTime1k_KV
    trj3d(t).fg1f = nan(3,length(vUseFrameIdx{t,1}));
    trj3d(t).fg2f = nan(3,length(vUseFrameIdx{t,1}));
    trj3d(t).hdf  = nan(3,length(vUseFrameIdx{t,1}));
    trj3d(t).jsTf = nan(3,length(vUseFrameIdx{t,1}));
    trj3d(t).jsBf = nan(3,length(vUseFrameIdx{t,1}));
    
    trj3d(t).fg1s = nan(3,length(vUseFrameIdx{t,1}));
    trj3d(t).fg2s = nan(3,length(vUseFrameIdx{t,1}));
    trj3d(t).hds  = nan(3,length(vUseFrameIdx{t,1}));
    trj3d(t).jsTs = nan(3,length(vUseFrameIdx{t,1}));
    trj3d(t).jsBs = nan(3,length(vUseFrameIdx{t,1}));
    
    if isstruct(fVideoInfo{t}) && isstruct(sVideoInfo{t})
        [fPath,fName,~] = fileparts(fVideoInfo{t}.path);
        [sPath,sName,~] = fileparts(sVideoInfo{t}.path);
        fIdx = cellfun(@(a) contains(a,fName), fList);
        sIdx = cellfun(@(a) contains(a,sName), sList);
        
        trj3d(t).fV = fullfile(fPath,fName); % front video path
        trj3d(t).sV = fullfile(sPath,sName); % side video path
        trj3d(t).vFrameT = vFrameTime{t}; % front/side video frame time points
        trj3d(t).vUseFrameIdx = vUseFrameIdx{t};
         
        
        if sum(fIdx)==1 && sum(sIdx)==1 && length(fron(fIdx).fg1IntX)==length(side(sIdx).fg1IntX) % in case, there are both front and side video files for this trial
            
            fFg1XY = [fron(fIdx).fg1IntX, fron(fIdx).fg1IntY]'; % finger1 front XY
            fFg2XY = [fron(fIdx).fg2IntX, fron(fIdx).fg2IntY]'; % finger2 front XY
            fHdXY  = [fron(fIdx).hdIntX, fron(fIdx).hdIntY]'; % hand front XY
            fjsTXY = [fron(fIdx).jsTIntX, fron(fIdx).jsTIntY]'; % joystick1 front XY
            fjsBXY = [fron(fIdx).jsBIntX, fron(fIdx).jsBIntY]'; % joystick2 front XY
            trj3d(t).jsfstPtFB = max([fron(fIdx).jsBfstPtX, fron(fIdx).jsBfstPtY]); 
            trj3d(t).jsfstPtFT = max([fron(fIdx).jsTfstPtX, fron(fIdx).jsTfstPtY]); 
            trj3d(t).jsfstPtSB = max([side(sIdx).jsBfstPtX, side(sIdx).jsBfstPtY]); 
            trj3d(t).jsfstPtST = max([side(sIdx).jsTfstPtX, side(sIdx).jsTfstPtY]); 
            
            sFg1XY = [side(sIdx).fg1IntX, side(sIdx).fg1IntY]'; % finger1 side XY
            sFg2XY = [side(sIdx).fg2IntX, side(sIdx).fg2IntY]'; % finger2 side XY
            sHdXY  = [side(sIdx).hdIntX, side(sIdx).hdIntY]'; % hand side XY
            sjsTXY = [side(sIdx).jsTIntX, side(sIdx).jsTIntY]'; % joystick1 side XY
            sjsBXY = [side(sIdx).jsBIntX, side(sIdx).jsBIntY]'; % joystick2 side XY
            
            % stereo triangulation
            if sum(sum(isnan(fFg1XY)==true))==0 && sum(sum(isnan(sFg1XY)==true))==0 
                [trj3d(t).fg1f,trj3d(t).fg1s] = stereoTriangulation(fFg1XY, sFg1XY, calibParams);
            end
            
            if sum(sum(isnan(fFg2XY)==true))==0 && sum(sum(isnan(sFg2XY)==true))==0
                [trj3d(t).fg2f,trj3d(t).fg2s] = stereoTriangulation(fFg2XY, sFg2XY, calibParams);
            end
            
            if sum(sum(isnan(fHdXY)==true))==0 && sum(sum(isnan(sHdXY)==true))==0
                [trj3d(t).hdf,trj3d(t).hds] = stereoTriangulation(fHdXY, sHdXY, calibParams);
            end
            
            if sum(sum(isnan(fjsTXY)==true))==0 && sum(sum(isnan(sjsTXY)==true))==0
                [trj3d(t).jsTf,trj3d(t).jsTs] = stereoTriangulation(fjsTXY, sjsTXY, calibParams);
            end
            
            if sum(sum(isnan(fjsBXY)==true))==0 && sum(sum(isnan(sjsBXY)==true))==0
                [trj3d(t).jsBf,trj3d(t).jsBs] = stereoTriangulation(fjsBXY, sjsBXY, calibParams);
            end
            
            %figure; plot3(trj3d.fg1f(1,:),trj3d.fg1f(2,:),trj3d.fg1f(3,:))
            %figure; plot3(trj3d.hdf(1,:),trj3d.hdf(2,:),trj3d.hdf(3,:))
            %figure; plot3(trj3d.jsTf(1,:),trj3d.jsTf(2,:),trj3d.jsTf(3,:))
            %figure; plot3(trj3d.jsBf(1,:),trj3d.jsBf(2,:),trj3d.jsBf(3,:))
            
            % For fingers and hand, take the median and smooth
            trj3d(t).allPartsF(:,:,1)=trj3d(t).fg1f;
            trj3d(t).allPartsF(:,:,2)=trj3d(t).fg2f;
            trj3d(t).allPartsF(:,:,3)=trj3d(t).hdf;
            
            trj3d(t).allPartsS(:,:,1)=trj3d(t).fg1s;
            trj3d(t).allPartsS(:,:,2)=trj3d(t).fg2s;
            trj3d(t).allPartsS(:,:,3)=trj3d(t).hds;
            
            trj3d(t).allPartsMedF = nanmedian(trj3d(t).allPartsF,3);
            trj3d(t).allPartsMedS = nanmedian(trj3d(t).allPartsS,3);
            
            % to use sg filter get sgfiltFramelen
            if length(trj3d(t).allPartsMedF) >= 101 % don't apply too aggressive smoothing here 101 is too much
                sgfiltFramelen = 101;
            elseif length(trj3d(t).allPartsMedF) < 101
                if mod(length(trj3d(t).allPartsMedF),2)==0
                    sgfiltFramelen = length(trj3d(t).allPartsMedF)-1; % the frame length for sg filter needs to be an odd number
                else
                    sgfiltFramelen = length(trj3d(t).allPartsMedF);
                end
            end
            trj3d(t).allPartsMedSgFron = sgolayfilt(trj3d(t).allPartsMedF',5,33)';%sgfiltFramelen)';
            trj3d(t).allPartsMedSgSide = sgolayfilt(trj3d(t).allPartsMedS',5,33)';%sgfiltFramelen)';        
            % plot3(trj3d(t).allPartsMedSgFron(1,:),trj3d(t).allPartsMedSgFron(2,:),trj3d(t).allPartsMedSgFron(3,:));
            % plot3(trj3d(t).allPartsMedSgSide(1,:),trj3d(t).allPartsMedSgSide(2,:),trj3d(t).allPartsMedSgSide(3,:));
            
            % For joystick trajectories, just smooth using sgolayfilt
            trj3d(t).jsTSgFron = sgolayfilt(trj3d(t).jsTf',5,55)'; % for js trj just smooth
            trj3d(t).jsTSgSide = sgolayfilt(trj3d(t).jsTs',5,55)'; % for js trj just smooth
            
            trj3d(t).jsBSgFron = sgolayfilt(trj3d(t).jsBf',5,55)';
            trj3d(t).jsBSgSide = sgolayfilt(trj3d(t).jsBs',5,55)';
            
            % select the valid portion of each trajectory when a file contained multiple trials
            if length(vUseFrameIdx{t})>length(vFrameTime{t})
                trj3d(t).allPartsMedSgFron = trj3d(t).allPartsMedSgFron(:,vUseFrameIdx{t});
                trj3d(t).allPartsMedSgSide = trj3d(t).allPartsMedSgSide(:,vUseFrameIdx{t});
            end
            % arrange rows to get intuitive XYZ trajectories
            trj3d(t).allPartsMedSgFronXYZ(1,:) = trj3d(t).allPartsMedSgFron(1,:); % X to be AP hand movement
            trj3d(t).allPartsMedSgFronXYZ(2,:) = trj3d(t).allPartsMedSgFron(3,:); % Y to be ML hand movement
            trj3d(t).allPartsMedSgFronXYZ(3,:) = -trj3d(t).allPartsMedSgFron(2,:); % Z to be UpOrDown hand movement 
            
            trj3d(t).allPartsMedSgSideXYZ(1,:) = trj3d(t).allPartsMedSgSide(1,:); 
            trj3d(t).allPartsMedSgSideXYZ(2,:) = trj3d(t).allPartsMedSgSide(3,:); 
            trj3d(t).allPartsMedSgSideXYZ(3,:) = -trj3d(t).allPartsMedSgSide(2,:);     
        
            trj3d(t).jsTSgFronXYZ(1,:) = trj3d(t).jsTSgFron(1,vUseFrameIdx{t}); 
            trj3d(t).jsTSgFronXYZ(2,:) = trj3d(t).jsTSgFron(3,vUseFrameIdx{t}); 
            trj3d(t).jsTSgFronXYZ(3,:) = -trj3d(t).jsTSgFron(2,vUseFrameIdx{t}); 
            
            trj3d(t).jsTSgSideXYZ(1,:) = trj3d(t).jsTSgSide(1,vUseFrameIdx{t});
            trj3d(t).jsTSgSideXYZ(2,:) = trj3d(t).jsTSgSide(3,vUseFrameIdx{t});
            trj3d(t).jsTSgSideXYZ(3,:) = -trj3d(t).jsTSgSide(2,vUseFrameIdx{t});
            
            trj3d(t).jsBSgFronXYZ(1,:) = trj3d(t).jsBSgFron(1,vUseFrameIdx{t}); 
            trj3d(t).jsBSgFronXYZ(2,:) = trj3d(t).jsBSgFron(3,vUseFrameIdx{t}); 
            trj3d(t).jsBSgFronXYZ(3,:) = -trj3d(t).jsBSgFron(2,vUseFrameIdx{t}); 
            
            trj3d(t).jsBSgSideXYZ(1,:) = trj3d(t).jsBSgSide(1,vUseFrameIdx{t});
            trj3d(t).jsBSgSideXYZ(2,:) = trj3d(t).jsBSgSide(3,vUseFrameIdx{t});
            trj3d(t).jsBSgSideXYZ(3,:) = -trj3d(t).jsBSgSide(2,vUseFrameIdx{t});
        end
    end
    clearvars tmp*
    fprintf('processed trial #%d\n', t);
end
clearvars t
save(fullfile(filePath,'trj3d.mat'), 'trj3d'); % save the raw pixel values

end

% rewardTrs = find([jsTime1k_KV(:).rewarded]==1);
% rwdTrI = 160;
% % x = trj3d(rewardTrs(rwdTrI)).allPartsMedSgSideXYZ(1,:); %trj3d(rewardTrs(rwdTrI)).vUseFrameIdx);
% % y = trj3d(rewardTrs(rwdTrI)).allPartsMedSgSideXYZ(2,:); %trj3d(rewardTrs(rwdTrI)).vUseFrameIdx);
% % z = trj3d(rewardTrs(rwdTrI)).allPartsMedSgSideXYZ(3,:); %trj3d(rewardTrs(rwdTrI)).vUseFrameIdx);
% % c = 1:sum(trj3d(rewardTrs(rwdTrI)).vUseFrameIdx); % generate a colormap;
% % 
% % figure;
% % patch([x nan],[y nan],[z nan],[c nan],'FaceColor','none','EdgeColor','interp')
% % colormap parula
% % colorbar
% % caxis([400 1000])
% 
% figure; hold on
% xT = trj3d(rewardTrs(rwdTrI)).jsTSgSideXYZ(1,trj3d(rewardTrs(rwdTrI)).vUseFrameIdx); %trj3d(rewardTrs(rwdTrI)).vUseFrameIdx);
% yT = trj3d(rewardTrs(rwdTrI)).jsTSgSideXYZ(2,trj3d(rewardTrs(rwdTrI)).vUseFrameIdx); %trj3d(rewardTrs(rwdTrI)).vUseFrameIdx);
% zT = trj3d(rewardTrs(rwdTrI)).jsTSgSideXYZ(3,trj3d(rewardTrs(rwdTrI)).vUseFrameIdx); %trj3d(rewardTrs(rwdTrI)).vUseFrameIdx);
% c = 1:sum(trj3d(rewardTrs(rwdTrI)).vUseFrameIdx); % generate a colormap;
% 
% patch([xT nan],[yT nan],[zT nan],[c nan],'FaceColor','none','EdgeColor','interp')
% 
% xB = trj3d(rewardTrs(rwdTrI)).jsBSgSideXYZ(1,trj3d(rewardTrs(rwdTrI)).vUseFrameIdx); %trj3d(rewardTrs(rwdTrI)).vUseFrameIdx);
% yB = trj3d(rewardTrs(rwdTrI)).jsBSgSideXYZ(2,trj3d(rewardTrs(rwdTrI)).vUseFrameIdx); %trj3d(rewardTrs(rwdTrI)).vUseFrameIdx);
% zB = trj3d(rewardTrs(rwdTrI)).jsBSgSideXYZ(3,trj3d(rewardTrs(rwdTrI)).vUseFrameIdx); %trj3d(rewardTrs(rwdTrI)).vUseFrameIdx);
% c = 1:sum(trj3d(rewardTrs(rwdTrI)).vUseFrameIdx); % generate a colormap;
% 
% patch([xB nan],[yB nan],[zB nan],[c nan],'FaceColor','none','EdgeColor','interp')
% colormap cool
% colorbar
% caxis([400 1000])
% 
% print(fullfile(filePath,'Figure',sprintf('tr#%d',rewardTrs(rwdTrI))),'-dpdf','-painters','-bestfit')
% jsTime1k_KV(rewardTrs(rwdTrI)).fVideo
% jsTime1k_KV(rewardTrs(rwdTrI)).sVideo

% fronFg1Fig2Dist = cell2mat(cellfun(@(a,b,c,d) sqrt((a-b).^2+(c-d).^2), {fron(:).fg1X}, {fron(:).fg2X}, {fron(:).fg1Y}, {fron(:).fg2Y},'Un',0)'); % front, point-by-point distance between two fingers
% fronFg1hdDist = cell2mat(cellfun(@(a,b,c,d) sqrt((a-b).^2+(c-d).^2), {fron(:).fg1X}, {fron(:).hdX}, {fron(:).fg1Y}, {fron(:).hdY},'Un',0)'); % front, point-by-point distance between two fingers

function [fTrj3d, sTrj3d] = stereoTriangulation(fXY, sXY, pr)
% performs triangulation using the stereo-triangulation data (pr) 

[fTrj3d, sTrj3d] = stereo_triangulation(fXY,sXY,pr.om,pr.T,pr.fc_left,pr.cc_left,pr.kc_left,...
pr.alpha_c_left,pr.fc_right,pr.cc_right,pr.kc_right,pr.alpha_c_right);

%plot3(fTrj3d(1,:),fTrj3d(2,:),fTrj3d(3,:))
%plot3(sTrj3d(1,:),sTrj3d(2,:),sTrj3d(3,:))

end


function  [intTrj,fstPt] = interplh(trj,trjlh,lm, zeroFirstNaNs)
% interpolates a timeseries (trj) based on it's likelihood (trjlh), enforce the Trjs to be within the frame dimension (lm)
%trj = tmpTrj.fg1Y;

%trj = tmpTrj.jsTX; trjlh = tmpTrj.jsTlh; lm = wd; zeroFirstNaNs = true;
intTrj = nan(length(trj),1);
x = 1:length(trj);

trj(trjlh<.9)=nan;
fstPt = NaN;

if zeroFirstNaNs
    isnanThres = .9; % in case for Joystick, there's a lot of NaNs expected while Js positioning
else
    isnanThres = .4; 
end

if sum(isnan(trj))/length(trj)<isnanThres % interpolation would be problematic with NaNs at the end
    valPts = strfind(num2str(trjlh>.9)','111');
    if sum(isnan(trj(end-4:end)))==0 % deal with last NaNs
        lastValPt = length(trj);
    elseif ~isempty(valPts)
        lastValPt = valPts(end)+2; % last valid point
    end
   
    fstPt = find(trjlh>.9,1,'first');
    tempTrj = trj(fstPt:lastValPt); % the valid portion of the trajectory
    tempTrjIdx = 1:length(trj)>=fstPt & 1:length(trj)<=lastValPt; % index for valid points
    if sum(isnan(tempTrj))>0
        % intrapolate points between the first and last valid points
        tempTrj(isnan(tempTrj)) = interp1(x(~isnan(tempTrj)),tempTrj(~isnan(tempTrj)),x(isnan(tempTrj)),'pchip');
    end
    trj(tempTrjIdx) = tempTrj;
    % extrapolate if each of the invalid portions at both sides is less than .1 of the full trj length
    if fstPt < length(trj)*.1 && lastValPt > length(trj)*.9
        trj(tempTrjIdx==0) = interp1(x(tempTrjIdx),trj(tempTrjIdx),x(tempTrjIdx==0),'pchip','extrap');
    end
    % put zeros for first NaNs, which often is the case for joystick trajectories
    if zeroFirstNaNs
        trj(1:fstPt) = 0;
    end
    trj = min(trj,lm);
    trj = max(trj,0);
    intTrj = trj;
end

end


% function  [intTrj,fstPt] = interplh(trj,trjlh,lm, dealFirstNaNs)
% % interpolates a timeseries (trj) based on it's likelihood (trjlh), enforce the Trjs to be within the frame dimension (lm)  
% %trj = tmpTrj.fg1Y;
% 
% %trj = tmpTrj.jsTX; trj1h = tmpTrj.jsTlh; lm = wd; dealFirstNaNs = true;
% intTrj = nan(length(trj),1); 
% x = 1:length(trj);
% 
% trj(trjlh<.9)=nan;
% fstPt = NaN;

%trjSpl = trj; 
%trjSpl(trjlh<.1)=nan;

%trjPch = trj; 
%trjPch(trjlh<.1)=nan;

%sum(isnan(trj))/length(trj); 

% if sum(isnan(trj(end-4:end)))==0 && sum(isnan(trj))/length(trj)<.4 % interpolation would be problematic with NaNs at the end
%     
%     if dealFirstNaNs % initial NaNs cannot be interpolated (e.g. joystick before/during positioning) - just interpolate subsequent to the 1st non-NaN value
%         fstPt = find(trjlh>.9,1,'first'); 
%         tempTrj = trj(fstPt:end); 
%         if sum(isnan(tempTrj))>0
%             x = 1:length(tempTrj); 
%             tempTrj(isnan(tempTrj)) = interp1(x(~isnan(tempTrj)),tempTrj(~isnan(tempTrj)),x(isnan(tempTrj)),'pchip');
%         end
%         
%         trj(fstPt:end) = tempTrj; 
%         trj = min(trj,lm); 
%         trj = max(trj,0); 
%         trj(1:fstPt) = 0;
%         intTrj = trj; 
%     %trjwnan = trj;
%     %figure; hold on;
%     %trj(isnan(trj)) = interp1(x(~isnan(trj)),trj(~isnan(trj)),x(isnan(trj)),'spline');
%     %trjSpl(isnan(trj)) = interp1(x(~isnan(trj)),trj(~isnan(trj)),x(isnan(trj)),'spline');
%     %trjPch(isnan(trj)) = interp1(x(~isnan(trj)),trj(~isnan(trj)),x(isnan(trj)),'pchip');
%     else % initial NaNs are not really an issue for body parts
%         trj(isnan(trj)) = interp1(x(~isnan(trj)),trj(~isnan(trj)),x(isnan(trj)),'pchip');
%         if sum(trj>lm)+sum(trj<0)==0 % enforce the trajectories within the frame
%             intTrj = trj;
%         end    
%         fstPt = NaN;
%     %hold on; 
%     %plot(trj); plot(trjSpl); plot(trjPch); 
%     %plot(trj)
%     %plot(trjwnan)
%     %hold off;
%     end
%     
% 
% end
% 
% end
