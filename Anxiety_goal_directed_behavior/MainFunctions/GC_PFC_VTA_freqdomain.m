%Granger-Causality analysis using the MVGC toolbox

clear all; clear functions; clc

%%
addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab\Functions'));
addpath(genpath('C:\Documents and Settings\jup36\Desktop\Matlab'));

%% load signal and timestamps modify this!
filedirectory = 'C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT\LFP'; cd(filedirectory); % directory for evt timestamps, behavioral data (tstpbeh)
load('LFP_AGB_periISPK_SpkLFP_PFCVTA_wolfpcutoff_102715','SpkLFP','tstpbeh'); % load LFP timeseries and behavioral time stamps      

filedirectory1 = 'C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT'; cd(filedirectory1); % directory for evt timestamps, behavioral data (tstpbeh)
load('PFC_sua_dat_090315','dat'); % load LFP timeseries and behavioral time stamps   

%% set initial parameters
uv.binsize = 0.001;    % binsize in sec (1 msec, point process)
uv.win = [2.25 2.25];  % peri-event time window (-1 to 0 s with each event occurring at time = 0, 100 ms extra-period added to both ends for binning)
uv.edges = [-2.25:uv.binsize:2.249];    % 1 ms bins to be used for histc
uv.detrend_win = [0.5, 0.05];    %detrending window size and step size 
uv.Fs = 1000;          %sampling frequency 

%% protocols for MVGC toolbox
gc.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
gc.icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

gc.morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
gc.momax     = 30;     % maximum model order for model order estimation

gc.acmaxlags = [];     % maximum autocovariance lags (empty for automatic calculation)
gc.acdectol = [];      % autocovariance decay tolerance
gc.regmode = [];

gc.tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
gc.alpha     = 0.05;   % significance level for significance test
gc.mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

gc.fs        = 1000;   % sample rate (Hz)
gc.fres      = [];     % frequency resolution (empty for automatic calculation)

gc.seed      = 0;      % random seed (0 for unseeded)

%% Multivariate Granger Causality analysis
PVfdgc = cell(15,4);   % a cell for the frequency domain granger causality: id, block 1, 2, 3 
VPfdgc = cell(15,4);   % a cell for the frequency domain granger causality: id, block 1, 2, 3

permPVfdgc = cell(15,4);   % a cell for the frequency domain granger causality: id, block 1, 2, 3 
permVPfdgc = cell(15,4);   % a cell for the frequency domain granger causality: id, block 1, 2, 3

counter = 0; % Counter
for f = 1:length(SpkLFP) % # of files
    
    if ~isempty(SpkLFP(f).Vispklfp) && ~isempty(SpkLFP(f).Pispklfp) % if there's valid corresponding PFC and VTA LFPs in both regions
        counter = counter + 1;
        % Detrend PFC & VTA lfp
        Plfpdt = locdetrend(SpkLFP(f).Pispklfp(1:4500,:), uv.Fs, uv.detrend_win); % detrend lfp; signal, sampling frequency, moving window
        Vlfpdt = locdetrend(SpkLFP(f).Vispklfp(1:4500,:), uv.Fs, uv.detrend_win); % detrend lfp; signal, sampling frequency, moving window
                
        %t = (1:size(Plfpdt, 1))'/uv.Fs; % time axis
        
        if strcmp(SpkLFP(f).PFCid(11),'A')==1     % Ascending order block
            tempidx = [1:size(SpkLFP(f).Pispklfp,2)];
            B1idx = tempidx<=50;           % B1 index
            B2idx = ~B1idx & tempidx<=100; % B2 index
            B3idx = ~B1idx & ~B2idx;       % B3 index
            
        elseif strcmp(SpkLFP(f).PFCid(11),'D')==1 % Descending order block
            tempidx = [1:size(SpkLFP(f).Pispklfp,2)];
            B3idx = tempidx<=50;           % B3 index
            B2idx = ~B3idx & tempidx<=100; % B2 index
            B1idx = ~B2idx & ~B3idx;       % B1 index
        end
        
        for b = 1:3 % block 1,2,3
            
            switch b
                case 1
                    idx = B1idx;
                    PVfdgc{counter,b} = SpkLFP(f).PFCid;    % put the file id in the 1st column PFC to VTA
                    VPfdgc{counter,b} = SpkLFP(f).VTAid;    % put the file id in the 1st column VTA to PFC
                    permPVfdgc{counter,b} = SpkLFP(f).PFCid; % 95% confidence interval bound
                    permVPfdgc{counter,b} = SpkLFP(f).VTAid; % 95% confidence interval bound
                case 2
                    idx = B2idx;
                case 3
                    idx = B3idx;
            end
            
            X = nan(2,length(Plfpdt),sum(idx)); % Region x Timepoints x Trials
            X(1,:,:) = reshape(Plfpdt(:,idx),[1,length(Plfpdt),sum(idx)]); % PFC lfp
            X(2,:,:) = reshape(Vlfpdt(:,idx),[1,length(Vlfpdt),sum(idx)]); % VTA lfp       
            
            gc.ntrials   = size(X,3);   % number of trials
            gc.nobs      = size(X,2);   % number of observations per trial (time in ms)
            
            % Estimate the model order
            [AIC,~,moAIC,~] = tsdata_to_infocrit(X,gc.momax,gc.icregmode);
            morder = moAIC;
            
            % VAR model estimation. Estimate VAR model of selected order from data.
            [A,SIG] = tsdata_to_var(X,morder,gc.regmode);   % A:VAR coefficients matrix, SIG:residuals covariance matrix
            
            % Check for failed regression
            assert(~isbad(A),'VAR estimation failed');
            
            [G,info] = var_to_autocov(A,SIG,gc.acmaxlags); % G:autocovariance sequence
            var_info(info,true); % report results (and bail out on error)
            
            %         F = autocov_to_pwcgc(G);
            %         assert(~isbad(F,false),'GC calculation failed');
            %
            %         % Significance test using theoretical null distribution, adjusting for multiple
            %         % hypotheses.
            %         pval = mvgc_pval(F,morder,gc.nobs,gc.ntrials,1,1,size(X,1)-2,gc.tstat); % take careful note of arguments!
            %         sig  = significance(pval,gc.alpha,gc.mhtc);
            
            % Frequency domain GC
            fd = autocov_to_spwcgc(G,500);
            assert(~isbad(fd,false),'spectral GC calculation failed');
            PVfdgc{counter,b+1} = squeeze(fd(2,1,:)); % put the frequency domain GC values in the 2nd column
            VPfdgc{counter,b+1} = squeeze(fd(1,2,:)); % put the frequency domain GC values in the 2nd column
            
            % Frequency domain GC empirical distribution by permutation
            fP = permtest_tsdata_to_spwcgc(X,morder,500,morder,1000); % 1st dim in fP: samples, 2nd dim: target (causee) variable, 3rd dim: the source (causal) variable, 4th dim: the fourth frequency.
            tempPVfdgc = sortrows(squeeze(fP(:,2,1,:)),-1); % permutation empirical distribution PFC to VTA
            tempVPfdgc = sortrows(squeeze(fP(:,1,2,:)),-1); % permutation empirical distribution VTA to PFC      
            
            permPVfdgc{counter,b+1} = tempPVfdgc(1:50,:); % 95% confidence interval bound
            permVPfdgc{counter,b+1} = tempVPfdgc(1:50,:); % 95% confidence interval bound
            
        end
        clearvars temp*
    else
        % Nothing happens for animals missing the VTA pair 
    end

end

%% Save data
cd('C:\Users\jup36\Desktop\Data\Anxiety_goal_directed_behavior\MAT\LFP\GC')
savename = 'GC_PFC_VTA_freqdomain_121815';
save(savename,'PVfdgc','VPfdgc','permPVfdgc','permVPfdgc')

%% plot 
ipPVfdgc = cell(size(PVfdgc,1),size(PVfdgc,2)); % preallocate the cell to put interpolated GC values
ipPVfdgc(:,1) = PVfdgc(:,1); % use the same file ID

xip = [1:501/501:501]; %[1:500/2000:500];  % the interpolation frequency space 


for f = 1:size(PVfdgc,1) % # of files
    for b = 2:4 % # of blocks
        int = 500/length(PVfdgc{f,b}); % original interval, which varies across all cases, 500 is the Nyquist frequency given the sampling rate 1000 Hz
        x = zeros(1,length(PVfdgc{f,b})); % preallcate x
        for i = 1:length(PVfdgc{f,b})  % # of data points
            x(i) = int*i; % frequency space using the current interval            
        end
        ipPVfdgc{f,b} = interp1(x, PVfdgc{f,b}, xip);
    end
end


ipVPfdgc = cell(size(VPfdgc,1),size(VPfdgc,2)); % preallocate the cell to put interpolated GC values
ipVPfdgc(:,1) = VPfdgc(:,1); % use the same file ID

xip = [1:501/501:501];  % the interpolation frequency space

for f = 1:size(VPfdgc,1) % # of files
    for b = 2:4 % # of blocks
        int = 500/length(VPfdgc{f,b}); % original interval, which varies across all cases, 500 is the Nyquist frequency given the sampling rate 1000 Hz
        x = zeros(1,length(VPfdgc{f,b})); % preallcate x
        for i = 1:length(VPfdgc{f,b})  % # of data points
            x(i) = int*i; % frequency space using the current interval            
        end
        ipVPfdgc{f,b} = interp1(x, VPfdgc{f,b}, xip);
    end
end


%% plot
cmap(2,:)=[0 114 189]./255; % pfc blue 
cmap(1,:)=[217 83 25]./255; % vta orange

xaxis = (xip>=1 & xip<=50);

% Block 1 
B1VPfdgcallfreq = cell2mat(ipVPfdgc(1:end,2)); % cell2mat of the all frequency band GCs
B1VPfdgc = nan(length(ipVPfdgc),25); % frequency bin adjusted average across four frequency bins
for b = 1:26
    B1VPfdgc(:,b) = nanmean(B1VPfdgcallfreq(:,2*(b-1)+1:2*b),2); % average across two freq bins 
end

B1PVfdgcallfreq = cell2mat(ipPVfdgc(1:end,2)); % cell2mat of the all frequency band GCs
B1PVfdgc = nan(length(ipPVfdgc),25); % frequency bin adjusted average across four frequency bins
for b = 1:26
    B1PVfdgc(:,b) = nanmean(B1PVfdgcallfreq(:,2*(b-1)+1:2*b),2); % average across two freq bins 
end

[mB1VPgc,~,sB1VPgc] = meanstdsem(cell2mat(ipVPfdgc(1:end,2))); % original g-causality mean across animals
[mB1PVgc,~,sB1PVgc] = meanstdsem(cell2mat(ipPVfdgc(1:end,2))); % original g-causality mean across animals

permBndPVb1 = zeros(length(permPVfdgc),length(permPVfdgc{1,2})); % conf interval at alpha = .001 
permBndVPb1 = zeros(length(permVPfdgc),length(permPVfdgc{1,2})); % conf interval at alpha = .001

permBndPVb2 = zeros(length(permPVfdgc),length(permPVfdgc{1,2})); % conf interval at alpha = .001 
permBndVPb2 = zeros(length(permVPfdgc),length(permPVfdgc{1,2})); % conf interval at alpha = .001

permBndPVb3 = zeros(length(permPVfdgc),length(permPVfdgc{1,2})); % conf interval at alpha = .001 
permBndVPb3 = zeros(length(permVPfdgc),length(permPVfdgc{1,2})); % conf interval at alpha = .001

% collect .001 bound from the permutation g-causality cell
for f = 1:length(permPVfdgc)
    permBndPVb1(f,:) = permPVfdgc{f,2}(1,:); % collect b1 .001 bound
    permBndVPb1(f,:) = permVPfdgc{f,2}(1,:); % collect b1 .001 bound
    
    permBndPVb2(f,:) = permPVfdgc{f,3}(1,:); % collect b1 .001 bound
    permBndVPb2(f,:) = permVPfdgc{f,3}(1,:); % collect b1 .001 bound
    
    permBndPVb3(f,:) = permPVfdgc{f,4}(1,:); % collect b1 .001 bound
    permBndVPb3(f,:) = permVPfdgc{f,4}(1,:); % collect b1 .001 bound
end

[permB1VPgc,~,permsB1VPgc] = meanstdsem(permBndVPb1); % permutation .1 % bound  
[permB1PVgc,~,permsB1PVgc] = meanstdsem(permBndPVb1); % permutation .1 % bound

figure(1);
boundedline(xip(xaxis), mB1VPgc(1,xaxis), sB1VPgc(1,xaxis), xip(xaxis), mB1PVgc(1,xaxis), sB1PVgc(1,xaxis), 'alpha','cmap', cmap) 
hold on;
plot(xip(xaxis), permB1VPgc(1,xaxis), 'Color', [217 83 25]./255, 'LineWidth', 2) 
plot(xip(xaxis), permB1PVgc(1,xaxis), 'Color', [0 114 189]./255, 'LineWidth', 2) 
xlim([1 50])
ylim([0 0.2])
set(gca,'XTick',[0:10:50],'FontSize',14);
set(gca,'YTick',[0:0.05:0.2],'FontSize',14);
set(gca,'TickDir','out')
hold off;

% Block 2
B2VPfdgcallfreq = cell2mat(ipVPfdgc(1:end,3)); % cell2mat of the all frequency band GCs
B2VPfdgc = nan(length(ipVPfdgc),25); % frequency bin adjusted average across four frequency bins
for b = 1:26
    B2VPfdgc(:,b) = nanmean(B2VPfdgcallfreq(:,2*(b-1)+1:2*b),2); % average across two freq bins 
end

B2PVfdgcallfreq = cell2mat(ipPVfdgc(1:end,3)); % cell2mat of the all frequency band GCs
B2PVfdgc = nan(length(ipPVfdgc),25); % frequency bin adjusted average across four frequency bins
for b = 1:26
    B2PVfdgc(:,b) = nanmean(B2PVfdgcallfreq(:,2*(b-1)+1:2*b),2); % average across two freq bins 
end

[mB2VPgc,~,sB2VPgc] = meanstdsem(cell2mat(ipVPfdgc(1:end,3)));
[mB2PVgc,~,sB2PVgc] = meanstdsem(cell2mat(ipPVfdgc(1:end,3)));

[permB2VPgc,~,permsB2VPgc] = meanstdsem(permBndVPb2); % permutation .1 % bound  
[permB2PVgc,~,permsB2PVgc] = meanstdsem(permBndPVb2); % permutation .1 % bound

figure(2);
boundedline(xip(xaxis), mB2VPgc(1,xaxis), sB2VPgc(1,xaxis), xip(xaxis), mB2PVgc(1,xaxis), sB2PVgc(1,xaxis), 'alpha','cmap', cmap) 
hold on;
plot(xip(xaxis), permB2VPgc(1,xaxis), 'Color', [217 83 25]./255, 'LineWidth', 2) 
plot(xip(xaxis), permB2PVgc(1,xaxis), 'Color', [0 114 189]./255, 'LineWidth', 2) 
xlim([1 50])
ylim([0 0.2])
set(gca,'XTick',[0:10:50],'FontSize',14);
set(gca,'YTick',[0:0.05:0.2],'FontSize',14);
set(gca,'TickDir','out')
hold off;

% Block 3
B3VPfdgcallfreq = cell2mat(ipVPfdgc(1:end,4)); % cell2mat of the all frequency band GCs
B3VPfdgc = nan(length(ipVPfdgc),25); % frequency bin adjusted average across four frequency bins
for b = 1:26
    B3VPfdgc(:,b) = nanmean(B3VPfdgcallfreq(:,2*(b-1)+1:2*b),2); % average across two freq bins 
end

B3PVfdgcallfreq = cell2mat(ipPVfdgc(1:end,4)); % cell2mat of the all frequency band GCs
B3PVfdgc = nan(length(ipPVfdgc),25); % frequency bin adjusted average across four frequency bins
for b = 1:26
    B3PVfdgc(:,b) = nanmean(B3PVfdgcallfreq(:,2*(b-1)+1:2*b),2); % average across two freq bins 
end

[mB3VPgc,~,sB3VPgc] = meanstdsem(cell2mat(ipVPfdgc(2:end,4)));
[mB3PVgc,~,sB3PVgc] = meanstdsem(cell2mat(ipPVfdgc(2:end,4)));

[permB3VPgc,~,permsB3VPgc] = meanstdsem(permBndVPb3); % permutation .1 % bound  
[permB3PVgc,~,permsB3PVgc] = meanstdsem(permBndPVb3); % permutation .1 % bound

figure(3);
boundedline(xip(xaxis), mB3VPgc(1,xaxis), sB3VPgc(1,xaxis), xip(xaxis), mB3PVgc(1,xaxis), sB3PVgc(1,xaxis), 'alpha','cmap', cmap) 
hold on;
plot(xip(xaxis), permB3VPgc(1,xaxis), 'Color', [217 83 25]./255, 'LineWidth', 2) 
plot(xip(xaxis), permB3PVgc(1,xaxis), 'Color', [0 114 189]./255, 'LineWidth', 2) 
xlim([1 50])
ylim([0 0.2])
set(gca,'XTick',[0:10:50],'FontSize',14);
set(gca,'YTick',[0:0.05:0.2],'FontSize',14);
set(gca,'TickDir','out')
hold off;





