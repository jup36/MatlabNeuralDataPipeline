function [PrevData] = TNC_SSPL_NevPreviewer(fileNameStr)

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

plotLogic =1;
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
            PrevData.unit(count).ts  = round(double(data.Data.Spikes.Timestamps(validInds) ./ 30));
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

if plotLogic
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
        title(['e' num2str(PrevData.unit(k).el) '-u' num2str(PrevData.unit(k).un)])
    end
end
%% CHECK STABILITY OF THE UNIT OVER TIME WITH THE SDF

% if plotLogic
%     figure(201); clf;
%     set(gcf,'Color','w');
%     for m=1:count
% 
%         numStamps                   = numel(PrevData.unit(m).ts);
%         delta                       = zeros(1,ceil(PrevData.unit(m).ts(numStamps)));
%         delta(PrevData.unit(m).ts)  = 1;
%         tmpSmooth                   = conv(delta,currParams.filter.kernel,'same') .* 1000;
% 
%         plot(tmpSmooth+((m-1).*10),'Color',cMap(m,:)); hold on;
%         axis off;
% 
%     end
% end
