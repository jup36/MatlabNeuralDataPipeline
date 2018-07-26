function [] = TNC_SSPL_Ns5Previewer(fileNameStr,chunk,numShanks,numSegsToView,arrayType)

%% GENERAL >>>  INITIALIZE PARAMETERS
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 100;
currParams.filter.causal       = 1;

[currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

if currParams.filter.causal
    currParams.filter.kernel(1:currParams.smthParams.decay.*15) = 0;
end
currParams.offset = 0;

currParams.winParams.prior     = 2e3;
currParams.winParams.after     = 5e3; 

plotLogic = 0;
PopData.currParams = currParams;

disp(' ');
disp('________________________________');
disp(' ');
disp('Initialized analysis parameters.');
disp('________________________________');
disp(' ');

%% DETERMINE THE NUMBER OF CHUNKS

    Ns5DATA = openNSx('report',fileNameStr);
    numSegs = floor(Ns5DATA.MetaTags.Duration ./ chunk);
    samples = floor(numSegs ./ numSegsToView);
    disp(' ');
    disp('________________________________________________________________');
    disp(' ');
    disp(['Data will be loaded and processed as ' num2str(numSegs) ' x ' num2str(chunk) ' seconds long segments.']);
    disp('________________________________________________________________');
    disp(' ');

%% EXAMINE THE EVENTS ON EACH CHANNEL (USE WAVELET FILTERING FOR FAST PERFORMANCE)
    
for shank = 1:numShanks

    for k = samples:samples:numSegs

        % LOAD ALL CHANNEL DATA
            timeStr = ['t:' num2str((k-1)*chunk) ':' num2str(k*chunk)];
                disp(' ');
                disp('________________________________________________________________');
                disp(' ');
                disp(['Loading data over the range of ' timeStr]);
                disp('________________________________________________________________');
                disp(' ');
            Ns5DATA = openNSx('report','read',fileNameStr,timeStr,'sec');

        % SMOOTH/FILTER DATA
            avgSeg      = median(Ns5DATA.Data,1);
            numChan     = size(Ns5DATA.Data,1); 
            cMap        = TNC_CreateRBColormap(8,'bo');

        % SOME CHECKING OF DATA FOR PROBES WHERE MULTIPLE SITES ARE USED FOR SORTING
            if ~strcmp(arrayType,'SingleSites')
                if numChan<64 & numChan>32
                    for zz=1:numChan
                        allRecChan(zz) = str2double(Ns5DATA.MetaTags.ElecLabel(zz,5:8));
                    end
                else
                    allRecChan = 1:numChan;
                end
            else
                allRecChan = 1:numChan;
            end
            
            if min(allRecChan) > 64
                allRecChan=allRecChan-64;
            end
            
        % GET THE ELECTRODE MAPPINGS
            [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(1,arrayType);


    %     for i=tmpChan:tmpChan
        for i = [1:8]
            
            thisChan = find(allRecChan==electrodeMatrix(i,shank));
            
            if numel( thisChan ) > 0
            
                % FILTER DATA FOR DISPLAY AND QUALITY CHECK
                rawData = sgolayfilt(Ns5DATA.Data(thisChan,:)-avgSeg,9,21);        

                % USE A WAVELET FILTERING METHOD FOR SPEED
                wname = 'db4'; 
                maxlevel = 5; % 7.5kHz low pass
                [c,l] = wavedec(rawData, maxlevel, wname);
                c = wthcoef('a', c, l);
                filtData = waverec(c, l, wname);

                % DISPLAY SOME 'RAW' DATA
                if plotLogic
                    figure(400); clf;
                    plot(Ns5DATA.Data(i,:),'b'); hold on;
                    plot(rawData,'k');
                    plot(filtData,'r');
                else
                    figure(400+shank); subplot(1,numSegsToView,k./samples);
                    if i==1
                        hold off;
                    end
                    plot(decimate(filtData,3)+(i.*750),'Color',cMap(i,:)); hold on;
                    axis([0 numel(decimate(filtData,3)) 0 6500]); axis off;
                    title(['sh' num2str(shank) ' - ' num2str(((k/2)*chunk)./60) ' min']);
                end
                
            else
                
                disp(['***ALERT*** Channel ' num2str(electrodeMatrix(i,shank)) ' was not recorded in this segment...']);
                
            end

        end

    end  
end        
        
        
        
        
        