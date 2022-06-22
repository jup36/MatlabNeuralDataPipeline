
%% PARAMETERS OF ANALYSIS
% Create the filename list from the current directory
% fileDirectory = 'C:\Users\babita\Documents\HHMI\Mitopark\Open Field\MP1\'
% fileDirectory = '/Users/dudmanj/Documents/Work/Janelia/_PROJECTS/MOVER_CODE/OPEN_FIELD_ANALYSIS/'
fileDirectory = '/Users/dudmanj/Documents/Work/Janelia/_PROJECTS/APPROACHINGNATURAL/'
% fileDirectory = './'
% Number of samples in running median
runMedWin   = 5;
 
%% DATA EXTRACTION
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(['________________Running RetrieveAllBehaviorData... Script_________________________'])
disp(' ');
disp(['Directory: ' fileDirectory]);
disp(['Timestamp: ' datestr(now)]);
disp(' ');
 
% create a list of all the files
allFiles = dir(sprintf('%s*.csv',fileDirectory));
numFiles = size(allFiles,1);

for k = 1:numFiles
 
    fileNameBase = allFiles(k).name(1:strfind(allFiles(k).name,'.')-1);
    dateStart = strfind(fileNameBase,'201');
    data.year = fileNameBase(dateStart(1):dateStart(1)+3);
    data.month = fileNameBase(dateStart(1)+4:dateStart(1)+5);
    data.day = fileNameBase(dateStart(1)+6:dateStart(1)+7);
    if data.month>12
        data.month = fileNameBase(dateStart(1)+4);
        data.day = fileNameBase(dateStart(1)+5:dateStart(1)+6);
    end
    date = datenum([data.month '-' data.day '-' data.year]);
    
%     data.age=age(fileNameBase,date);    
 
    % Script to load csv behavior data
    clear disp* trial* 
 
    skipInitTrials = 0;
 
    fileNameStr = [fileDirectory fileNameBase '.csv'];
    disp([fileDirectory fileNameBase '.csv']);
 
    tmpData = dlmread(fileNameStr,',',2,0);
    numSamps = size(tmpData,1);
    
    %___________________________________________________________
    % TIMESTAMPS (MS)
    data.tStamps = tmpData(:,1);
    
    %___________________________________________________________
    % SMOOTH AND NICE UP THE ACCEL DATA
    data.accel.x = tmpData(:,2);
    data.accel.y = tmpData(:,3);
    data.accel.z = tmpData(:,4);
    data.accel.mag=sqrt(data.accel.x.^2+data.accel.y.^2+data.accel.z.^2);
    
    %___________________________________________________________
    % SMOOTH AND NICE UP THE GYRO DATA
    data.gyro.x = tmpData(:,5);
    data.gyro.y = tmpData(:,6);
    data.gyro.z = tmpData(:,7);
    data.gyro.mag = sqrt(data.gyro.x.^2+data.gyro.y.^2+data.gyro.z.^2);
 
    %___________________________________________________________
    % DIGITAL SIGNALS FROM BEHAVIOR SYSTEM
    data.dig    = tmpData(:,8);
 
    %___________________________________________________________
    % SMOOTH AND NICE UP THE POSITION DATA
    data.posRaw.xA = tmpData(:,9);
    data.posRaw.yA = tmpData(:,10);
    data.posRaw.xB = tmpData(:,11);
    data.posRaw.yB = tmpData(:,12);
 
    %___________________________________________________________
    % POSITION DATA APPEARS TO BE OVERSAMPLED AND SO I NEED TO REMOVE ZEROS
    % AND THEN CREATE THE CONTINUOUS FUNCTION AGAIN WITH SMOOTHING
    data.posSmth.xA = sgolayfilt(data.posRaw.xA,11,31);
    data.posSmth.yA = sgolayfilt(data.posRaw.yA,11,31);
    data.posSmth.xB = sgolayfilt(data.posRaw.xB,11,31);
    data.posSmth.yB = sgolayfilt(data.posRaw.yB,11,31);
 
    data.posSmth.x  = (data.posSmth.xA+data.posSmth.xB)./2;
    data.posSmth.y  = (data.posSmth.yA+data.posSmth.yB)./2;    
    
    %___________________________________________________________
    % FURTHER SMOOTH POSITION DATA BY CREATING RUNNING MEDIANS    
    data.posSmth.Xrmed=medfilt1(data.posSmth.x,runMedWin);
    data.posSmth.Yrmed=medfilt1(data.posSmth.y,runMedWin);
            
    %___________________________________________________________
    % EUCLIDEAN DISTANCE WITH SMOOTHED POSITION DATA
    i=2:numel(data.posSmth.Xrmed);
    data.velocity.vel(i) = sqrt( (data.posSmth.Xrmed(i)-data.posSmth.Xrmed(i-1)).^2 + (data.posSmth.Yrmed(i)-data.posSmth.Yrmed(i-1)).^2 );
 
    %___________________________________________________________
    % FURTHER SMOOTH VELOCITY DATA 
    data.velocity.velRM     = abs(sgolayfilt(data.velocity.vel,5,11));
    data.velocity.velRML(i) = real(log(data.velocity.velRM(i)));
    
   %___________________________________________________________
   % TAKE CUMULATIVE SUMS OF DISTANCE TO GET DISTANCE FROM ORIGIN
    data.velocity.cumdistance = (cumsum(data.velocity.vel))';
    
        % Write the data structure to a specific name
    disp(' ');
    disp(' ');
    disp([fileNameBase '.data = data;']);
    eval([fileNameBase '.data = data;']);
    clear data;
    
    disp(['Output data in the structure: ' fileNameBase]);
    disp(' ');
   
    clear tmpData;
 
end


disp(['________________Completed Script_______________________________________________________'])
disp(' ');
disp(' ');
disp(' ');
disp(' ');
