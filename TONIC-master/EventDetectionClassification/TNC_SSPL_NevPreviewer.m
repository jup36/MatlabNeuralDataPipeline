function [PrevData] = TNC_SSPL_NevPreviewer(fileNameStr,alignChan)

%% GENERAL >>>  INITIALIZE PARAMETERS
currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 10;
currParams.filter.causal       = 1;

[currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);

if currParams.filter.causal
    currParams.filter.kernel(1:currParams.smthParams.decay.*15) = 0;
end
currParams.offset = 0;

currParams.winParams.prior     = 1e2;
currParams.winParams.after     = 1e3; 

plotLogic   =1;
plotLogicC  =0;
rasterLogic =1;

PopData.currParams = currParams;

disp(' ');
disp('________________________________');
disp(' ');
disp('Initialized analysis parameters.');
disp('________________________________');
disp(' ');

%% LOAD THE DATA
clear PrevData;
[data] = TNC_LoadData(0, 0, fileNameStr);

%% COLLECT ONLINE SORTED UNITS

onSort = find(data.Data.Spikes.Unit>0 & data.Data.Spikes.Electrode<128);

totPossUnits = max(unique(data.Data.Spikes.Unit));

onSortChan = unique(data.Data.Spikes.Electrode(onSort));

count = 1;

for i=1:numel(onSortChan)
    
    for j=1:totPossUnits
        
        validInds = find(data.Data.Spikes.Electrode == onSortChan(i) & data.Data.Spikes.Unit==j);
        if numel(validInds)>2
            PrevData.unit(count).wfs = double(data.Data.Spikes.Waveform(validInds,:));
            PrevData.unit(count).ts  = ceil(double(data.Data.Spikes.Timestamps(validInds) ./ 30));
            PrevData.unit(count).el  = double(onSortChan(i));
            PrevData.unit(count).un  = double(j);
            count = count+1;
        end
        
    end
    

end

count = count-1;
if count<1
    disp(' ');
    disp('________________________________');
    disp(' ');
    disp('There are no online sorted units.');
    disp('________________________________');
    disp(' ');

    PrevData = [];
    plotLogic = 0;
end

%% DISPLAY ONLINE SORTED UNITS

cMap = TNC_CreateRBColormap(count,'bo');

if plotLogic & count>1
    disp(' ');
    disp('________________________________');
    disp(' ');
    disp(['Plotting ' num2str(count) ' sorted units for ' fileNameStr]);
    disp('________________________________');
    disp(' ');

    cMap = TNC_CreateRBColormap(count,'bo');
    figure(200); clf;
    set(gcf,'Color','w');
    for k=1:count
        subplot(ceil(count./8),8,k);
        shadedErrorBar(1:numel(PrevData.unit(k).wfs(1,:)) , mean(PrevData.unit(k).wfs) , std(PrevData.unit(k).wfs,[],1) , {'-k','color',cMap(k,:)});
        axis([0 50 -1800 1800]); axis off;
        title([num2str(k) '>> e' num2str(PrevData.unit(k).el) '-u' num2str(PrevData.unit(k).un)])
    end
end

%% DISPLAY RASTERS ALIGNED TO CHOSEN CHANNEL
figure(199); clf;
set(gcf,'Color','w');

if rasterLogic & count>1 & alignChan~=0

    alignInds = find(data.Data.Spikes.Electrode==137 & data.Data.Spikes.Unit==1);
    times = round(double(data.Data.Spikes.Timestamps(alignInds))./30);
    
    for m=1:count
        
        valids      = find(PrevData.unit(m).ts<max(times)+currParams.winParams.after+1000);
        delta       = zeros(1,max(times)+currParams.winParams.after+1000);
        delta(PrevData.unit(m).ts(valids)+1) = 1;
        tmpSmooth   = conv(delta,currParams.filter.kernel,'same');
        tmpSmoothZ  = ( tmpSmooth - mean(tmpSmooth) ) ./ std( tmpSmooth );

        [alignData] = TNC_AlignRasters(tmpSmoothZ , PrevData.unit(m).ts , -1 ,times, [currParams.winParams.prior currParams.winParams.after],0,1);
        PrevData.unit(m).csAligned    = alignData.image.psthAVG;        
        PrevData.unit(m).csAlignedE   = alignData.image.psthSEM;        

        subplot(ceil(count./8),8,m);
        shadedErrorBar([-currParams.winParams.prior:currParams.winParams.after] , PrevData.unit(m).csAligned , PrevData.unit(m).csAlignedE , {'-k','color',cMap(m,:)});
        axis([-currParams.winParams.prior currParams.winParams.after min(PrevData.unit(m).csAligned)*1.25 max(PrevData.unit(m).csAligned)*1.25]); axis off;
        title(['e' num2str(PrevData.unit(m).el) '-u' num2str(PrevData.unit(m).un)])

        if m==2
            fid = fopen('yokedata.csv','w');
            for p=1:numel(times)
                yokeInds = find(PrevData.unit(m).ts>times(p)+40 & PrevData.unit(m).ts<times(p)+200);
                forYoke.trial(p).ts = round(PrevData.unit(m).ts(yokeInds) - times(p));
                forYoke.numSpks(p) = numel(yokeInds);
                fprintf(fid,'%d,',p);
                fprintf(fid,'%d,',numel(yokeInds));
                for r=1:numel(yokeInds)
                    if r==numel(yokeInds)
                        fprintf(fid,'%d',forYoke.trial(p).ts(r));
                    else
                        fprintf(fid,'%d,',forYoke.trial(p).ts(r));
                    end
                end
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
        
    end
    
end

%% CHECK STABILITY OF THE UNIT OVER TIME WITH THE SDF

if plotLogicC & count>1
    figure(201); clf;
    set(gcf,'Color','w');
    for m=1:count

        numStamps                   = numel(PrevData.unit(m).ts);
        delta                       = zeros(1,ceil(PrevData.unit(m).ts(numStamps)));
        delta(PrevData.unit(m).ts)  = 1;
        tmpSmooth                   = conv(delta,currParams.filter.kernel,'same') .* 1000;

        plot(tmpSmooth+((m-1).*10),'Color',cMap(m,:)); hold on;
        axis off;

    end
end

