%% Trajectory analysis code
% Goals:
% Illustrate nice joystick trajectories by transforming:
% 1. Separate out scaling (amplitude) from direction (angle)
% 2. Analyze main axis and "off axis"

%% Preprocess movement data

%=============================
% USE ELLIPTICAL FILTER
par.sr = 25000;
par.detect_fmin = 1;
par.detect_fmax = 1000;

[spk_map] = TNC_CreateRBColormap(8,'mbr');
spk_kernel = TNC_CreateGaussian(500,25,1000,1);
lk_kernel = TNC_CreateGaussian(500,60,1000,1);

lick_logic = zeros(1,numel(lick));
lick_shift = [0 lick(1:end-1)];
all_licks = find(lick>0.0075 & lick_shift<0.0075);
figure(2); clf; plot(lick); hold on; plot(all_licks,ones(1,numel(all_licks)).*0.0075,'r*');
lick_logic(all_licks) = 1;

figure(2); hold off;

position.x = TNC_FilterData2(Xpos,par);
position.y = TNC_FilterData2(Ypos,par);

% plot(Xpos); hold on;
plot(position.x);hold on;
plot(position.y);
% plot(reachNoStim.currEvt{3}(:,1),ones(numel(reward.currEvt{3}(:,1)),1)*-0.1,'ro');
 
tmp = find(position.x<-0.1 & position.y<-0.1);
tmp2 = tmp(find(diff(tmp)>1));
tmp3 = tmp2(1:numel(tmp2)-1);
plot(tmp3,-0.1*ones(1,numel(tmp3)),'bs');

Ypos_cln = position.y;
Xpos_cln = position.x;

for qq=1:numel(tmp3)
    val = Ypos_cln(tmp3(qq)-60);
    Ypos_cln([tmp3(qq)-60:tmp3(qq)+30]) = val*ones(1,91);
    val = Xpos_cln(tmp3(qq)-60);
    Xpos_cln([tmp3(qq)-60:tmp3(qq)+30]) = val*ones(1,91);
end
position.y_cln = sgolayfilt(Ypos_cln,3,201);
position.x_cln = sgolayfilt(Xpos_cln,3,201);
[tmp_th,position.r] = cart2pol(position.x_cln,position.y_cln);
[tmp_th,position.r_v] = cart2pol([0 diff(position.x_cln)],[0 diff(position.y_cln)]);
position.th = conv(tmp_th,spk_kernel,'same');

figure(3); clf;
plot(position.y_cln); hold on;
plot(position.y);
plot(reward.currEvt{3}(:,1),ones(numel(reward.currEvt{3}(:,1)),1)*-0.1,'ro');
% [rt_mvmt] = TNC_TransformTrajectory(position.x_cln(1,:),position.y_cln(1,:),50,24);

shift = 0;
position.d_cln = sqrt( position.x_cln.^2 + position.y_cln.^2 );
[sink.x] = TNC_ExtTrigWins(position.x_cln,rewardNoStim.currEvt{3}(:,1)+shift,[3000 2000]);
[sink.y] = TNC_ExtTrigWins(position.y_cln,rewardNoStim.currEvt{3}(:,1)+shift,[3000 2000]);
[sink.d] = TNC_ExtTrigWins(position.d_cln,rewardNoStim.currEvt{3}(:,1)+shift,[3000 2000]);
[sink.s] = TNC_ExtTrigWins(conv((1000.*[0 diff(position.d_cln)]),spk_kernel,'same'),rewardNoStim.currEvt{3}(:,1)+shift,[3000 2000]);
[sink.a] = TNC_ExtTrigWins([0 diff(sgolayfilt([0 diff(position.d_cln)],3,101))],rewardNoStim.currEvt{3}(:,1)+shift,[3000 2000]);
[sink.l] = TNC_ExtTrigWins(conv(lick_logic,lk_kernel,'same'),rewardNoStim.currEvt{3}(:,1)+shift,[3000 2000]);
[sink.th] = TNC_ExtTrigWins(position.th,rewardNoStim.currEvt{3}(:,1)+shift,[3000 2000]);
[sink.r] = TNC_ExtTrigWins(position.r,rewardNoStim.currEvt{3}(:,1)+shift,[3000 2000]);

num_trials = size(sink.x.wins,1);

%% Extract movements and rotate for kinematics extraction
lick_kern = TNC_CreateGaussian(500,50,1000,1);
position.d_cln = sqrt( position.x_cln.^2 + position.y_cln.^2 );
figure(10); clf;
plot(position.d_cln); hold on; plot(position.y_cln);
lick_df = conv(lick_logic,lick_kern,'same');

before = 1000;
after = 2000;
range = 1:rewardNoStim.currEvt{3}(end,1)+5000;
thr_range = 0.005:0.005:0.04; 

% Find all movements
for ev_thr = thr_range
    [events] = TNC_ExtractPeaks(position.d_cln(range),ev_thr,1000,1);
    
    [sink.r] = TNC_ExtTrigWins(abs(sole),events.starts,[before after]);
    clear traj;
    for kk=1:size(sink.r.wins,1)
        if max(sink.r.wins(kk,1500:2500))>0.1
            traj.rew(kk) = 1;
        else
            traj.rew(kk) = 0;
        end
    end    
    chk_thr.rew(find(ev_thr==thr_range)) = sum(traj.rew)./size(rewardNoStim.currEvt{3},1);
    
end

figure(100); clf; plot(thr_range,chk_thr.rew,'ko'); axis([min(thr_range) max(thr_range) 0 1]); box off; grid on;

[~,best_thr] = min(abs(chk_thr.rew-1));
ev_thr = thr_range(best_thr);

[events] = TNC_ExtractPeaks(position.d_cln(range),ev_thr,500,1);

% Align data to movements
[sink.x] = TNC_ExtTrigWins(position.x_cln,events.starts,[before after]);
[sink.y] = TNC_ExtTrigWins(position.y_cln,events.starts,[before after]);
[sink.l] = TNC_ExtTrigWins(lick_df,events.starts,[before after]);
[sink.r] = TNC_ExtTrigWins(abs(sole),events.starts,[before after]);
[sink.la] = TNC_ExtTrigWins(laser,events.starts,[before after]);

figure(11); clf;
subplot(151); imagesc(sink.x.wins);
subplot(152); imagesc(sink.y.wins);
subplot(153); imagesc(sink.r.wins);
subplot(154); imagesc(sink.l.wins);
subplot(155); imagesc(sink.la.wins);

% extract joystick trajectories and do the transform to get them in the
% same orientation / isolate scaling.

[~,maxi] = max(abs(sink.y.avg));
range = before-50:maxi;
% range = before-50:2000;

figure(12); clf;
method = 'hull'; clear traj

for kk=1:size(sink.x.wins,1)

    switch method
        case 'proc'
            [d,Z,transform] = procrustes([sink.x.avg(range)' sink.y.avg(range)'] , [sink.x.wins(kk,range)' sink.y.wins(kk,range)']);
            scale(kk) = transform.b;
            rot_traj = [sink.x.wins(kk,range)' sink.y.wins(kk,range)']*transform.T;
        case 'hull'
            [theta,radius] = cart2pol(sink.x.wins(kk,range),sink.y.wins(kk,range));
            theta = theta-pi/2;
            [scale(kk),maxdi] = max(radius);            
            rot_traj = [sink.x.wins(kk,range)' sink.y.wins(kk,range)']*[cos(theta(maxdi)) -sin(theta(maxdi)) ; sin(theta(maxdi)) cos(theta(maxdi))] ;
    end
    traj.x(kk,:) = rot_traj(:,1)'-rot_traj(1,1);
    traj.y(kk,:) = rot_traj(:,2)'-rot_traj(1,2);
    traj.th(kk) = theta(maxdi);
    traj.amp(kk) = sqrt( (rot_traj(maxdi,1)-rot_traj(1,1)).^2 + (rot_traj(maxdi,2)-rot_traj(1,2)).^2 );
    
    if max(sink.r.wins(kk,1500:2500))>0.1
        traj.rew(kk) = 1;
    else
        traj.rew(kk) = 0;
    end
    if max(sink.la.wins(kk,500:2000))>0.1
        traj.laze(kk) = 1;
    else
        traj.laze(kk) = 0;
    end

    subplot(211);
%     plot(sink.x.wins(kk,range)+(0.1*kk),sink.y.wins(kk,range),'color',[0.5 0.5 0.5 0.9]); hold on; 
%     plot(sink.x.wins(kk,range(maxdi))+(0.1*kk),sink.y.wins(kk,range(maxdi)),'r.'); hold on; 
    plot(sink.x.wins(kk,range),sink.y.wins(kk,range),'color',[0.5 0.5 0.5 0.1]); hold on; 
    plot(sink.x.wins(kk,range(maxdi)),sink.y.wins(kk,range(maxdi)),'r.'); hold on; 
    subplot(212);
%     plot(rot_traj(:,1)+(0.1*kk),rot_traj(:,2)-rot_traj(1,2),'color',[0.85 0 0 0.3]); hold on; 
%     plot(rot_traj(maxdi,1)+(0.1*kk),rot_traj(maxdi,2)-rot_traj(1,2),'r.'); hold on; 
    plot(rot_traj(1:maxdi,1)-rot_traj(1,1),rot_traj(1:maxdi,2)-rot_traj(1,2),'color',[0.85 0 0 0.1]); hold on; 
    plot(rot_traj(maxdi,1)-rot_traj(1,1),rot_traj(maxdi,2)-rot_traj(1,2),'r.'); hold on; 
%     plot(rot_traj(:,1)+(0.1*kk),rot_traj(:,2),'color',[0 0.42 0.85 0.3]); hold on; 
%     plot(rot_traj(maxdi,1)+(0.1*kk),rot_traj(maxdi,2),'.','color',[0 0.42 0.85]); hold on; 
    
end
subplot(211);
plot(sink.x.avg(range),sink.y.avg(range),'k','linewidth',2);
axis([-0.1 0.1 -0.1 0.1]); box off; grid on;
subplot(212);
% plot([0 14],[0 0],'k-');
% plot(mean(traj.x),mean(traj.y),'k','linewidth',2);
axis([-0.1 0.1 -0.1 0.1]); box off; grid on;

disp(['Found ' num2str(sum(traj.rew)) ' rewards.']);

%% Compute the SDF 
clear physDat;
spk_kern = TNC_CreateGaussian(500,25,1000,1);

for qq=1:size(spkTimesCellCTX,2)
    
    times = round(double(spkTimesCellCTX{1,qq}));
    physDat.str_logic(qq) = spkTimesCellCTX{5,qq};
    physDat.depth(qq) = spkTimesCellCTX{4,qq}(1,2);
    
    delta = zeros(1,numel(laser));
    delta(times(times<numel(laser) & times>0)) = 1;
    
    physDat.sdf(qq,:) = conv(delta,spk_kern,'same');
    physDat.sdfz(qq,:) = (physDat.sdf(qq,:) - mean(physDat.sdf(qq,:))) ./ std(physDat.sdf(qq,:));
    
    fprintf('processed cell # %d\n', qq) % report unit progression
end

[mA, m] = compute_mapping(physDat.sdfz', 'PCA', 10); % m.M has the PC loadings (mapping)

[sink_pop] = TNC_ExtTrigWins3d(physDat.sdfz,events.starts(traj.rew==1 & traj.laze==0  & traj.amp>ev_thr),[before after]);
[sink_pcs] = TNC_ExtTrigWins3d(mA',events.starts(traj.rew==1 & traj.laze==0 & traj.amp>ev_thr),[before after]);

[~,sp_mi] = max(sink_pop.win_avg,[],2);
[~,sp_mi_s] = sort(sp_mi);

figure(13); clf;
imagesc(sink_pop.win_avg(sp_mi_s,:));

figure(14); clf;
plot(sink.y.range,sink_pcs.win_avg(1:2,:));

[deeps,deepi] = sort(physDat.depth);

figure(15); clf;
subplot(121);
imagesc(m.M(deepi,1:2));
subplot(122);
plot(physDat.depth(deepi),1:size(physDat.sdfz,1));
set(gca,'YDir','rev');
hold on;
plot([500 500],[0 250],'k');
plot([1000 1000],[0 250],'k');
% plot(zeros(1,numel(find(stimE.tagE(deepi)==-1))),find(stimE.tagE(deepi)==-1),'ro');

figure(11); clf;
subplot(151); imagesc([sink.x.wins(traj.laze==0 & traj.amp>ev_thr,:) ; sink.x.wins(traj.laze==1,:)]);
subplot(152); imagesc([sink.y.wins(traj.laze==0 & traj.amp>ev_thr,:) ; sink.y.wins(traj.laze==1,:)]);
subplot(153); imagesc([sink.r.wins(traj.laze==0 & traj.amp>ev_thr,:) ; sink.r.wins(traj.laze==1,:)]);
subplot(154); imagesc([sink.l.wins(traj.laze==0 & traj.amp>ev_thr,:) ; sink.l.wins(traj.laze==1,:)]);
subplot(155); imagesc([sink.la.wins(traj.laze==0 & traj.amp>ev_thr,:) ; sink.la.wins(traj.laze==1,:)]);

figure(20); clf;
perturb_dat = [traj.amp(traj.laze==0 & traj.amp>ev_thr) traj.amp(traj.laze==1 & traj.amp>ev_thr)];
perturb_grp = [zeros(1,numel(traj.amp(traj.laze==0 & traj.amp>ev_thr))) ones(1,numel(traj.amp(traj.laze==1 & traj.amp>ev_thr)))];
boxplot(perturb_dat,perturb_grp);

%% Solve for KN dimension
% Ideas:
% 1. dPCA to clear up movement timing component from the decision component
% 2. Use cleaned movement quantification to try and solve for KN
% 3. Use multiple regression (like decoding approach) to get cleaner KN dim (separated from theta dim)

[ord_map] = TNC_CreateRBColormap(6,'cpb');
[pm_map] = TNC_CreateRBColormap(12,'mbr');
[ex_map] = TNC_CreateRBColormap(8,'exag');
[sink_pop] = TNC_ExtTrigWins3d(physDat.sdfz,events.starts,[before after]);

% KN dim - solve beta values for 50 randperms of events
batch = 50;
for rr=1:50
    % get BATCH of trials
    valids = find(traj.laze==0);
    these_inds = valids(randperm(numel(valids),batch));
    
    % compute z-scored activity modulation for each trial
    cnt =1; clear mags;
    for ss=these_inds
       mags(:,cnt) = trapz(sink_pop.wins(:,before-100:before+400,ss),2);
       cnt=cnt+1;
    end
    spds = traj.amp(these_inds);
    
    % compute multiple regression
    betas(:,rr) = pinv(mags')*spds';
        
end

% average betas to get consensus KN dimension
figure(50); clf;
% imagesc([betas(sp_mi_s,:) zeros(250,1) mean(betas(sp_mi_s,:),2)]); colormap(ex_map);
kn.dim = mean(betas,2);
subplot(121);
imagesc(kn.dim(deepi),[-5e-6 5e-6]);
colormap(pm_map);
subplot(122);
plot(physDat.depth(deepi),1:size(physDat.sdfz,1));
set(gca,'YDir','rev');
hold on;
plot([500 500],[0 250],'k');
plot([1000 1000],[0 250],'k');


figure(29); clf;
h = histogram(traj.amp,6);

% First look at pop activity binned into quintiles of movement amplitude
figure(30); clf;
for pp=1:10
    [pc(pp).sink] = TNC_ExtTrigWins(mA(:,pp)',events.starts,[before after]);
end

[kn.sink] = TNC_ExtTrigWins((physDat.sdfz'*kn.dim)',events.starts,[before after]);
cnt = 1;
for qq=1:numel(h.BinEdges)-1
    
    inds = find(traj.amp>h.BinEdges(qq) & traj.amp<=h.BinEdges(qq+1));
    
    for pp=1:5
        figure(30);
        bin_pc(qq).dat(pp,:) = mean(pc(pp).sink.wins(inds,:));
        subplot(1,6,pp); 
        plot(pc(pp).sink.range,bin_pc(qq).dat(pp,:),'color',ord_map(qq,:),'linewidth',2); hold on;
        axis([-1000 2000 -10 15]);
    end
    
    subplot(1,6,6);
    plot(kn.sink.range,mean(kn.sink.wins(inds,:)),'color',ord_map(qq,:),'linewidth',2); hold on;
    axis([-1000 2000 -5e-5 20e-5]);

        
    kn.avg_amp(cnt) = mean(traj.amp(inds));
    kn.avg_load(cnt) = trapz(mean(kn.sink.wins(inds,before-100:before+400)));
    
    figure(31);
    subplot(1,numel(h.BinEdges)-1,qq); 
    imagesc(mean(sink_pop.wins(sp_mi_s,:,inds),3),[0 10]); colormap(ex_map);

    cnt = cnt+1;
end

figure(32); scatter(kn.avg_amp,kn.avg_load,100,1:numel(kn.avg_amp),'filled'); colormap(ord_map); box off; grid on; ylabel('KN loading'); xlabel('Mvmt. amp. (a.u.)');

%% Project tagged neurons onto KN dimension

laze_events = TNC_ExtractPeaks(-laser,1,500,1);

sdfz_ctx = physDat.sdfz; % physDat.sdfz(physDat.str_logic==0,:);
sdf_ctx = physDat.sdf; %physDat.sdf(physDat.str_logic==0,:);
kn_ctx = (sdfz_ctx'*kn.dim)';
kn.dim_ctx = kn.dim; % kn.dim(physDat.str_logic==0);
deep_ctx = physDat.depth; % physDat.depth(physDat.str_logic==0);

% Loading of tagged neurons onto KN dimension
[kn.laze] = TNC_ExtTrigWins(kn_ctx,laze_events.starts,[before after]);
[kn.sink_ctx] = TNC_ExtTrigWins(kn_ctx,events.starts,[before after]);

figure(40); clf;
plot(kn.sink_ctx.range,kn.sink_ctx.avg); hold on; plot(kn.sink.range,kn.sink.avg,'k'); plot(kn.laze.range,kn.laze.avg);

% Effect of inactivation on KN dimension
[tag_sink] = TNC_ExtTrigWins3d(sdfz_ctx,laze_events.starts,[before after]);
figure(41); clf;
% imagesc(tag_sink.win_avg);
put_tag = []; tag_thr = -0.35;
for kk=1:size(tag_sink.wins,1)
    tmp = mean( tag_sink.wins(kk,:,:) , 3 ) ;
    if mean(tmp(before+250:before+500))-mean(tmp(1:500)) < tag_thr & mean(tmp(before:before+250))-mean(tmp(1:500)) < tag_thr
        plot(tmp-mean(tmp(1:500)),'color',[0 0.45 0.85]); hold on;
        put_tag = [put_tag kk];
    else
        plot(tmp-mean(tmp(1:500)),'color',[0.5 0.5 0.5 0.1]);
    end
end

figure(42); clf; 
hist(kn.dim_ctx);
hold on;
plot(kn.dim_ctx(put_tag),ones(1,numel(put_tag)),'r*');

figure(44); clf;
hist(deep_ctx);
hold on;
plot(deep_ctx(put_tag),ones(1,numel(put_tag)),'r*');


[~,kn_is] = sort(kn.dim_ctx,'descend');

physDat.put_tag = put_tag;
[sink_pop_ctx] = TNC_ExtTrigWins3d(sdfz_ctx,events.starts(traj.rew==1 & traj.laze==0  & traj.amp>ev_thr),[before after]);
figure(43); clf;
for zz=1:size(sdf_ctx,1)
    if numel(find(zz==put_tag))==1
        plot(mean(sink_pop_ctx.wins(zz,:,:),3),'color',[0 0.45 0.85]); hold on;
    elseif numel(find(zz==kn_is(1:numel(put_tag))))==1
        plot(mean(sink_pop_ctx.wins(zz,:,:),3),'color',[0.85 0 0]); hold on;
    else
        plot(mean(sink_pop_ctx.wins(zz,:,:),3),'color',[0.5 0.5 0.5 0.1]); hold on;
    end
end


