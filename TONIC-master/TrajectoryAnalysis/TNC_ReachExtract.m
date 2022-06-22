function [reach] = TNC_ReachExtract( pathName, fileName, varargin )
eric_true = 1;

%% Reach parameters

csvLogic = 0;

if nargin==2

    % Required fields in the params struct:
    min_space       = 50; % generous initial threshold
    ampl_density_thr= 0.75;
    duration_thr    = 50;
    distance_thr    = 25;
    plot_frac       = 0.1;

    % PARAMETERS FOR MERGING
    ang_diff_thr    = pi./24;
    chg_ang_thr     = pi;
    gap_thr         = 500; % gap in time prior to reach that means, trivially, it is a new reach

    chan.ids        = [137,138,139,142,143];
    chan.labels    = {'x' 'y' 'lick' 'rew' 'thr'};

elseif nargin==3

    % Required fields in the params struct:
    min_space       = 50; % generous initial threshold
    ampl_density_thr= 0.75;
    duration_thr    = 50;
    distance_thr    = 25;
    plot_frac       = 0.1;

    % PARAMETERS FOR MERGING
    ang_diff_thr    = pi./24;
    chg_ang_thr     = pi;
    gap_thr         = 500; % gap in time prior to reach that means, trivially, it is a new reach  
    
    chan = varargin{1};

else
    
    chan = varargin{1};
    params = varargin{2};

    % Required fields in the params struct:
    min_space       = params.min_space; % generous initial threshold
    ampl_density_thr= params.ampl_density_thr;
    duration_thr    = params.duration_thr;
    distance_thr    = params.distance_thr;
    plot_frac       = params.plot_frac;

    % PARAMETERS FOR MERGING
    ang_diff_thr    = params.ang_diff_thr;
    chg_ang_thr     = params.chg_ang_thr;
    gap_thr         = params.gap_thr; % gap in time prior to reach that means, trivially, it is a new reach  
    
end

%% DATA LOADING

    dataRate = fileName(numel(fileName)-2:numel(fileName));
    
    switch dataRate
    
        case 'ns4'    

            disp(['Loading data from an ' dataRate ' file...']);
            [cont_data] = openNSx('read','report',[pathName fileName]);

%             for i=1:numel(cont_data.ElectrodesInfo)
%                 
%                 eval([ 'chan.' chan.labels{find(chan.ids == cont_data.ElectrodesInfo(i).ElectrodeID )} ' = ' num2str(i) ';']);
%                 disp([ 'chan.' chan.labels{find(chan.ids == cont_data.ElectrodesInfo(i).ElectrodeID )}  ' = ' num2str(i) ';']);
%                 
%             end

            disp(' ');
            for p=1:numel(cont_data.ElectrodesInfo)
                cont_data.ElectrodesInfo(p).ElectrodeID
                ind = find( chan.ids == cont_data.ElectrodesInfo(p).ElectrodeID );
                if numel(ind)==1
                    disp(['chan.' chan.labels{ind} ' = ' num2str(p) ';']);
                    eval(['chan.' chan.labels{ind} ' = ' num2str(p) ';']);
                else
                    disp(['chan.' chan.labels{p} ' = ' num2str(0) ';']);
                    eval(['chan.' chan.labels{p} ' = ' num2str(0) ';']);                    
                end
            end
            
            pos.x =  decimate(double(cont_data.Data(chan.x,:)), 10 );
            pos.y =  decimate(double(cont_data.Data(chan.y,:)), 10 );        

            metaData = cont_data.MetaTags;
            
        case 'ns2'

            pos.x =  double(cont_data.Data(chan.x,:));
            pos.y =  double(cont_data.Data(chan.y,:));

            disp(['Loading data from an ' dataRate ' file...']);
            [cont_data] = openNSx('read',[pathName fileName]);

            disp(' ');
            for p=1:numel(cont_data.ElectrodesInfo)
                ind = find( chan.ids == cont_data.ElectrodesInfo(p).ElectrodeID );
                eval(['chan.' chan.labels{ind} ' = ' num2str(p) ';']);
                disp(['chan.' chan.labels{ind} ' = ' num2str(p) ';']);
            end
            
            metaData = cont_data.MetaTags;
            
        case 'csv'

            csvLogic = 1;
            gap_thr = round(gap_thr ./ 35);
            duration_thr = round(duration_thr ./ 35);
            [data_struct] = TNC_LoadCSV_Mover([pathName fileName]);
            pos.x = data_struct.pos.x';
            pos.y = data_struct.pos.y';
            metaData = data_struct.metaData;

    end
    
%% PRE-PROCESS RAW DATA
    
    pos.xs = sgolayfilt( pos.x , 9 , 21 );
    pos.ys = sgolayfilt( pos.y , 9 , 21 );

    if ~csvLogic
        [countsx,centersx] = hist(pos.xs,1000);
        pos.xs = pos.xs - centersx( countsx==max(countsx) );
        [countsy,centersy] = hist(pos.ys,1000);
        pos.ys = pos.ys - centersy( countsy==max(countsy) );
        disp(['Estimating an offset of (' num2str(centersx( countsx==max(countsx) )) ' , ' num2str(centersy( countsy==max(countsy) )) ')']);
    end
    
    pos.vs = [0 sqrt( diff(pos.xs).^2 + diff(pos.ys).^2 )];
    if csvLogic
        pos.v_win = sgolayfilt( pos.vs , 3 , 11 );
    else
        pos.v_win = conv(pos.vs,[zeros(1,5) ones(1,10) zeros(1,5)]./10,'same');
    end
%     pos.v_win = moveavg(pos.vs,10);
    [theta,pos.rho] = cart2pol(pos.xs,pos.ys);
    
    x_diff = [0 diff(pos.xs)];
    y_diff = [0 diff(pos.ys)];
    pos.angle = atan2(sgolayfilt(y_diff,3,101),sgolayfilt(x_diff,3,101));

    plot_range          = 1:round(numel(pos.xs)*plot_frac);
    
    disp('...complete');
    disp(' ');
   
%% use hist to find the mode and then create a fake distribution around the mode using 1:mode reflected

    disp('Finding suprathreshold velocities to begin reach segmentation...');

    numSamples = numel(pos.v_win);
    [counts,centers] = hist(pos.v_win,2000);
    hi_cnt  = max(counts); % find the mode; define threshold based upon the low s.d. estimate of baseline
    low_c   = centers( find(counts<(hi_cnt./5), 1 ) );
    mode_c  = centers( counts==hi_cnt );
    hi_c    = centers( find(counts>(hi_cnt.*0.37) , 1 , 'last' ) );
    threshold = centers( find(cumsum(counts)./sum(counts) > 0.85 , 1) ); % no more than 20% of session should be intrareach
    string_thresh = centers( find(cumsum(counts)./sum(counts) > 0.7 , 1) );
    if threshold < ( mean(pos.v_win) + 2.*std(pos.v_win) ) && csvLogic==0
        threshold = mean(pos.v_win) + 2.*std(pos.v_win);
        string_thresh = mean(pos.v_win) + 1.*std(pos.v_win);
    end
    
% get suprathreshold velocities and find boundaries
    c_intra_reach       = find( pos.v_win > threshold );
    pos.angle            = pos.angle .* ( pos.v_win > threshold );
    frac_intra_reach    = numel(c_intra_reach) ./ numel(pos.v_win);
    if csvLogic
        inter_reach_loc     = find([3 diff(c_intra_reach)] > 2);
    else
        inter_reach_loc     = find([21 diff(c_intra_reach)] > 20);
    end
    numReaches        = numel( inter_reach_loc );
    sb_starts           = [ c_intra_reach( inter_reach_loc(1:numReaches-1) ) ];
    sb_stops            = [ c_intra_reach( inter_reach_loc(2:numReaches) - 1 ) ];
    pos.angle           = pos.angle .* (pos.v_win > string_thresh);

    figure(20); clf;
    plot(pos.vs(plot_range)); hold on; plot(pos.v_win(plot_range)); plot(pos.angle(plot_range).*10); plot(sqrt(pos.rho(plot_range)));
    plot( sb_starts(sb_starts<max(plot_range)), ones(1,numel(sb_starts(sb_starts<max(plot_range)))).*threshold, 'g^');
    plot( sb_stops(sb_stops<max(plot_range)), ones(1,numel(sb_stops(sb_stops<max(plot_range)))).*threshold, 'r*');
    drawnow;    

    disp('...complete');
    disp(' ');
    
%% Use thresholding to find putative starts   
    disp('Finding valid reach starts...');
    
    debug_on = 0;
    
    sc_starts           = sb_starts;
    sc_stops            = sb_stops;
    total_vel           = zeros(1,numel(sc_starts));
    initial_angle       = zeros(1,numel(sc_starts));
    outward             = zeros(1,numel(sc_starts));
    diff_angle          = zeros(1,numel(sc_starts));
    init_position       = zeros(1,numel(sc_starts));
    
    for j=2:numel(sb_starts)-1
        fix_start = find( pos.v_win(sb_stops(j-1)+1:sb_starts(j)) < string_thresh , 1 , 'last' );
        if numel(fix_start)==0
            sc_starts(j) = sb_starts(j);
        else
            sc_starts(j) = sb_stops(j-1) + fix_start;
        end

        fix_stop = find( pos.v_win(sb_stops(j):sb_starts(j+1)) < string_thresh , 1 , 'first' );
        if numel(fix_stop)==0
            sc_stops(j) = sb_stops(j);
        else
            sc_stops(j) = sb_stops(j) + fix_stop;
        end
    end
    
    for j=1:numel(sc_starts)
        angles          = pos.angle(sc_starts(j):sc_stops(j));
        velos           = pos.v_win(sc_starts(j):sc_stops(j));
        rhos            = pos.rho(sc_starts(j):sc_stops(j));
        total_vel(j)    = sum(velos);
        weights         = velos ./ sum(velos);
        initial_angle(j)= sum(angles.*weights) ./ sum(weights); % weight the angle estimate by inst. velocity
        return_angle(j) = atan2(0-pos.ys(sc_starts(j)) , 0-pos.xs(sc_starts(j)));
        outward(j)      = pos.rho(sc_stops(j)) - pos.rho(sc_starts(j));
        init_position(j)= sqrt(pos.xs(sc_starts(j)).^2 + pos.ys(sc_starts(j)).^2);        
        diff_angle(j)   = angleDiff(initial_angle(j), return_angle(j));
    end
    
    figure(3); clf;
    plot(pos.v_win(plot_range)); hold on; plot(pos.angle(plot_range).*25); 
    plot( sc_starts(sc_starts<max(plot_range)), initial_angle(sc_starts<max(plot_range)).*25, 'g^');
    
    figure(4); clf;
        subplot(221); rose(initial_angle); hold on; rose(return_angle); 
        subplot(222); scatter(initial_angle,outward,16,log(total_vel)); 
        subplot(223); scatter(log(init_position),abs(diff_angle),16,outward);
        subplot(224); hist(abs(diff_angle),314);
    
    if debug_on
        suspicious_reaches = find( abs(diff_angle) < ang_diff_thr );
        for jj=1:numel(suspicious_reaches)
            j = suspicious_reaches(jj);
            figure(5); plot(pos.xs(sc_starts(j):sc_stops(j)),pos.ys(sc_starts(j):sc_stops(j)),'k');
            axis([-750 750 -750 750]); pause(0.2);
        end
    end
        
    disp('...complete');
    disp(' ');

%% Decide how to merge starts into stops
% Previous step gives me lots of reach segments which I now have to decide
% about whether to merge or not

    disp('Merging / validating reach segments into reaches...');
    
%     walk through all of the reaches and look at whether 
%         (0) is this a new beginning of movement - distant enough from prior stop
%         (1) the current reach segment is valid - centripetal and fast enough
%         (2) the current reach should be merged with the next reach - change in angle, time between
        
    seg_ind = 2;
    mrg_ind = 1;
    cnt     = 0;
    clear reach new_reach;
    reach.start(1)= sc_starts(1);                
    reach.stop(1) = sc_stops(1);
    reach.start_ecc(1) = pos.rho( sc_starts(1) );
    dec_criteria = [];            
    h = waitbar(0,'Please wait...');
    
    while seg_ind < numel(sc_starts)

        pre_gap = sc_starts(seg_ind)-reach.stop(mrg_ind);
        
        if pre_gap > gap_thr % trivially independent
            
            if abs(diff_angle(seg_ind)) > ang_diff_thr % looks outward (>8 degrees diff then return trajectory), its a keeper                
                mrg_ind             = mrg_ind+1;
                reach.start(mrg_ind)= sc_starts(seg_ind);                
                reach.stop(mrg_ind) = sc_stops(seg_ind);
                reach.start_ecc(mrg_ind) = pos.rho( sc_starts(seg_ind) );
            end
            
        else % compare with the prior reach segment to see if they should be joined
            further_out     = pos.rho( sc_stops(seg_ind) ) - (reach.start_ecc(mrg_ind) .* 1.1);
            change_dir      = angleDiff( initial_angle(seg_ind) , initial_angle(seg_ind-1) );
            
            if abs(diff_angle(seg_ind)) > ang_diff_thr && further_out>0 % it is still going outward
                if change_dir < chg_ang_thr % if its roughly in same direction, then merge 
                    reach.stop(mrg_ind) = sc_stops(seg_ind); % update current reach with the stop of this segment (i.e. merge)
                else % start a new reach if there was a dramatic change in direction
                    mrg_ind             = mrg_ind+1;
                    reach.start(mrg_ind)= sc_starts(seg_ind);                
                    reach.stop(mrg_ind) = sc_stops(seg_ind);                    
                    reach.start_ecc(mrg_ind) = pos.rho( sc_starts(seg_ind) );
                end                
            else
                % just ignore, must be the return component, and wait until a big enough gap to start over
            end            
        end

        seg_ind = seg_ind+1;
        
    end
    
    reach.num = mrg_ind;
    
    disp(['...found ' num2str(mrg_ind) ' valid reaches from ' num2str(numel(sc_starts)) ' candidate segments.']);
    disp('...filtering reaches based upon amplitude and duration ');

    for i=1:reach.num
        rhos   = pos.rho(reach.start(i):reach.stop(i));
        ampl(i) = max(rhos) - pos.rho(1);        
        dist(i) = pos.rho(reach.stop(i)) - pos.rho(reach.start(i));
        durs(i) = numel(rhos);
        waitbar(i / reach.num);
    end
    
    figure(20); subplot(311); hist(log10(ampl),100); title('Max amplitude (log)'); subplot(312); hist(log10(abs(dist)),100); title('End distance'); subplot(313); hist(log10(durs),100); title('Duration');
    figure(21); subplot(211); plot(log10(ampl),log10(dist),'k.'); xlabel('amplitude'); ylabel('distance'); subplot(212);  plot(log10(durs),log10(dist),'k.');
    %ask user to pick boundary??

    cnt = 0;
        
    for i=1:reach.num
        if durs(i) > duration_thr
            cnt = cnt+1;
            new_reach.start(cnt)    = reach.start(i);
            new_reach.stop(cnt)     = reach.stop(i);
            new_reach.traj(cnt).x   = pos.xs(reach.start(i):reach.stop(i));
            new_reach.traj(cnt).y   = pos.ys(reach.start(i):reach.stop(i));
            new_reach.traj(cnt).v   = pos.vs(reach.start(i):reach.stop(i));
        end
        waitbar(i / reach.num);
    end
    
    disp(['...final estimate from JOSH: ' num2str(cnt) ' valid reaches.']);
    disp(' ');
    
    close(h);
    
    clear reach;
    
    reach            = new_reach;
    reach.pos     = pos;

    figure(3); clf;
    plot_inds = find(new_reach.stop<max(plot_range));
    subplot(211);
    plot(pos.v_win(plot_range)); hold on; plot(sqrt(pos.rho(plot_range))); plot(new_reach.start(plot_inds),pos.v_win(new_reach.start(plot_inds)),'g*');
    plot(new_reach.stop(plot_inds),pos.v_win(new_reach.stop(plot_inds)),'r*');
    subplot(212);
    plot(pos.xs(plot_range)); hold on; plot(new_reach.start(plot_inds),pos.xs(new_reach.start(plot_inds)),'g*');
    plot(new_reach.stop(plot_inds),pos.xs(new_reach.stop(plot_inds)),'r*');
    
%% GET REACH STATS ON FINAL SET
    for j=1:numel(reach.start)
        
        angles                  = pos.angle(reach.start(j):reach.stop(j));
        velos                   = pos.v_win(reach.start(j):reach.stop(j));
        rhos                    = pos.rho(reach.start(j):reach.stop(j));
        
        reach.stats.maxV(j)     = max(velos);
        reach.stats.avgV(j)     = mean(velos);
        reach.stats.avgA(j)     = mean(diff(velos));
        reach.stats.totA(j)     = sum(diff(velos));
        reach.stats.dist(j)     = trapz(velos);
        reach.stats.ampl(j)     = pos.rho(reach.stop(j)) - pos.rho(reach.start(j));       
                weights         = velos ./ sum(velos);
        reach.stats.angle(j)    = sum(angles.*weights) ./ sum(weights); % weight the angle estimate by inst. velocity
        
    end
    
%     [reach.stats.lowd, reach.stats.mapping] = compute_mapping([reach.stats.maxV' reach.stats.avgV' reach.stats.avgA' reach.stats.totA' reach.stats.dist' reach.stats.ampl' reach.stats.angle'],'tSNE');
    [reach.stats.lowd, reach.stats.mapping] = compute_mapping([zscore(reach.stats.maxV)' zscore(reach.stats.avgV)' zscore(reach.stats.avgA)' zscore(reach.stats.totA)' zscore(reach.stats.dist)' zscore(reach.stats.ampl)' zscore(reach.stats.angle)'],'tSNE');

    
%% GET EVENT DATA

    if ~csvLogic
        disp('Getting behavioral event data...');

        if chan.rew>0
            rew_cont_data   = double(cont_data.Data(chan.rew,:));
            [rew_inds]      = TNC_FindEventsInContData(rew_cont_data,min_space);

            reach.rew_inds = round(rew_inds./10);
            reach.thr_inds = reach.rew_inds - 1000;
        else
            reach.rew_inds = [];
            reach.thr_inds = [];
        end
        
        reach.dataRate = dataRate;
        reach.metaData = metaData;

        disp(['...found: ' num2str(numel(rew_inds)) ' trials.']);
        disp(' ');
    end
    
%% PROCESS WITH ERIC'S ANALYSIS FOR COMPARISON    
if eric_true
    disp('Extracting reaches with Eric algorithm...');

    [ reachStart reachStop reach0 pos1 pos2] = getReachTimes( [ pos.x ; pos.y ] );
    
    reach.yttri.reachStart  = reachStart;
    reach.yttri.reachStop   = reachStop;
    reach.pos.v_ma          = reach0;
              
    figure(6); clf; 
    plot(pos.rho); hold on; plot(pos.v_win,'k');
    plot(reachStart, reach0(reachStart),'g*','MarkerSize',20)
    plot(reachStop, reach0(reachStop),'r*','MarkerSize',20)
    axis([0 max(plot_range) 0 max(reach0)]);
    plot(reach.start, reach0(reach.start),'go','MarkerSize',20)
    plot(reach.stop, reach0(reach.stop),'ro','MarkerSize',20)

    figure(7); clf;
    subplot(121); imagesc(pos1);
    subplot(122); imagesc(pos2);

    psthWin = [1e3,3e3];
    [pos_rho]   = TNC_ExtTrigWins(reach.pos.rho,reach.start(reach.start>psthWin(1)),psthWin);
    [pos_rho_y] = TNC_ExtTrigWins(reach.pos.rho,reach.yttri.reachStart(reach.yttri.reachStart>psthWin(1)),psthWin);

    figure(5); clf;
    shadedErrorBar(-psthWin(1):psthWin(2),mean(pos_rho.wins,1),std(pos_rho.wins,[],1)./sqrt(size(pos_rho.wins,1)),{'Color',[0 0.67 1]});
    hold on;
    shadedErrorBar(-psthWin(1):psthWin(2),mean(pos_rho_y.wins,1),std(pos_rho_y.wins,[],1)./sqrt(size(pos_rho_y.wins,1)),{'Color',[1 0.67 0]});   

    disp(['...final estimate from ERIC: ' num2str(numel(reachStart)) ' valid reaches.']);
    disp(' ');
end
      
%% Rotate and re-sample reaches to get into a common reference
% 
% numReaches = numel(dataStructure.data.reach.start);
% ymean = [];
% xmean = [];
% angles = [];
% sampleRate = [];
% 
% for i=1:numReaches
% 
%     startInd = find(dataStructure.data.trajectories.contXYsmooth(:,1) >= dataStructure.data.reach.start(i),1);
%     stopInd = find(dataStructure.data.trajectories.contXYsmooth(:,1) >= dataStructure.data.reach.stop(i),1);
%     currTraj = [dataStructure.data.trajectories.contXYsmooth(startInd:stopInd,2) , dataStructure.data.trajectories.contXYsmooth(startInd:stopInd,3)];
%     currTraj(:,1) = currTraj(:,1) - currTraj(1,1);
%     currTraj(:,2) = currTraj(:,2) - currTraj(1,2);
%     
%     try
%         k = convhull(currTraj(:,1),currTraj(:,2));
%     catch ME
%         k = 1:numel(currTraj(:,1));
%     end
% 
%     hullPntsX = currTraj(k,1);
%     hullPntsY = currTraj(k,2);
%     hullDists = sqrt( hullPntsX.^2 + hullPntsY.^2 );
%     maxDisplace = find(hullDists == max(hullDists),1);                      % index of max displacement
%     maxDispTotalTime = find(currTraj(:,1)==hullPntsX(maxDisplace) & currTraj(:,2)==hullPntsY(maxDisplace), 1, 'first');     % index of max displacement in the entire trajectory    
%     
%     majReachAngle = atan2(currTraj(maxDispTotalTime,2)-currTraj(1,2) , currTraj(maxDispTotalTime,1)-currTraj(1,1));
%     
%     % rotate to lie along the y axis
%     rotAngle = (pi./2) - majReachAngle;
%     rotationTransform = [cos(rotAngle) -sin(rotAngle) ; sin(rotAngle) cos(rotAngle)];
%     newTraj = currTraj*rotationTransform';
% 
%     % fit a polynomial to the trajectory data
%     yPF = polyfit(1:maxDispTotalTime,newTraj(1:maxDispTotalTime,2)',10);
%     xPF = polyfit(1:maxDispTotalTime,newTraj(1:maxDispTotalTime,1)',10);    
%     resample = maxDispTotalTime ./ 100;
% 
%     % build matrix for average of all fits :: raw
%     raw.reach(i).trajectory = newTraj(1:maxDispTotalTime,:);
%     raw.reach(i).majReachAngle = majReachAngle;
%     
%     % build matrix for average of all fits :: re-scaled    
%     ymean = [ymean ; polyval(yPF,[resample:resample:maxDispTotalTime]+1)];
%     xmean = [xmean ; polyval(xPF,[resample:resample:maxDispTotalTime]+1)];
%     angles = [angles ; majReachAngle];
%     sampleRate = [sampleRate ; resample];    
% 
% %     figure(figNum); 
% %     if i==1
% %         clf;
% %     end
% %     subplot(221);
% %     plot(currTraj(1:maxDispTotalTime,1),currTraj(1:maxDispTotalTime,2),'k',[-75 75],[0 0],'k--',[0 0],[-75 75],'k--'); axis tight; hold on;
% %     subplot(221);
% %     plot(newTraj(1:maxDispTotalTime,1),newTraj(1:maxDispTotalTime,2),'k',[-30 30],[0 0],'k--',[0 0],[-25 75],'k--'); axis tight; hold on;%     subplot(223);
% %     plot(1:100,polyval(yPF,resample:resample:maxDispTotalTime),'r',[0 100],[0 0],'k--'); axis tight; hold on;
% %     subplot(224);
% %     plot(1:100,polyval(xPF,resample:resample:maxDispTotalTime),'r',[0 100],[0 0],'k--'); axis tight; hold on;
% % 
% %     subplot(222);
% % %     plot(polyval(xPF,1:maxDispTotalTime),polyval(yPF,1:maxDispTotalTime),'r',[-30 30],[0 0],'k--',[0 0],[-25 75],'k--'); axis tight; hold on;
% %     plot(1:maxDispTotalTime,newTraj(1:maxDispTotalTime,2)); hold on; %,'r',[-30 30],[0 0],'k--',[0 0],[-25 75],'k--'); axis tight; hold on;
%     
% end
% 
% figure(figNum); clf;
% subplot(121);
% % shadedErrorBar(1:100,mean(xmean,1),std(xmean,[],1),{'k','Color',[1 0 0]}); hold on;
% plot(1:100,std(xmean,[],1),'r'); hold on;
% shadedErrorBar(1:100,mean(ymean,1),std(ymean,[],1)./sqrt(size(ymean,1)-1),{'k','Color',[0 0.55 1]});
% plot([0 100],[0 0],'k--'); axis([0 100 -5 25]); ylabel('Distance (a.u.)'); xlabel('Reach phase (%)');
% subplot(122);
% plot(mean(xmean,1),mean(ymean,1),'k','LineWidth',2); hold on;
% plot([-10 10],[0 0],'k--',[0 0],[-5 20],'k--'); axis tight;
% 
% dataStructure.xmean     = xmean;
% dataStructure.ymean     = ymean;
% dataStructure.angles    = angles;
% dataStructure.sampleRate= sampleRate;
% dataStructure.raw       = raw;

      