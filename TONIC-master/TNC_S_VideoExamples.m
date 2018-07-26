%% PRELIMINARIES

session.config.vid.ch = 144;


%% LOAD THE DATA

% video files
vidFile = 'FullBody-m38-20110324.seq';              % video file name
[seq_info, fid] = TNC_ReadSeqHeader(vidFile);       % TNC_ReadSeqHeader extracts information about the video
vidFile2 = 'Whiskers-m38-20110324-002.seq';         % vidoe file name
[seq_info2, fid2] = TNC_ReadSeqHeader(vidFile2);    % TNC_ReadSeqHeader extracts information about the video

% neural data (?)
base = 'm38-20110324-buz278b-005';
data = openNEV(['/' base '.nev'],'read','report','nosave'); % blackrock spike data (?)

%% GET FRAME STAMPS FROM VIDEO DATA
ts_frames = double( data.Data.Spikes.TimeStamp( find(data.Data.Spikes.Electrode==session.config.vid.ch) ) ); % getting frame timestamps
disp(['Total frames: ' num2str(seq_info.NumberFrames) ' ; ' num2str(seq_info2.NumberFrames) ' | Total stamps: '  num2str(numel(ts_frames))]);

%% PROCESS MOVEMENTS FROM CONTINUOUS DATA AND FIND REACH ONSETS

ns4_data = openNSx(['/' base '.ns4'],'read','report');

session.config.position.x.ch = 137;
session.config.position.y.ch = 138;
session.config.water_port.ch = 139;
session.config.lick.ch       = 139;
session.config.events.reward.ch = 142;
session.config.events.threshold.ch = 143;

chan.ids       = [session.config.position.x.ch, session.config.position.y.ch, session.config.water_port.ch, session.config.events.reward.ch, session.config.events.threshold.ch];
chan.labels    = {'x' 'y' 'lick' 'rew' 'thr'};

[reach] = TNC_ReachExtract( './' , [base '.ns4'], chan );

%% LOAD THE PHYSIOLOGY DATA TO PLOT ALIGNED PSTHS WITH THE VIDEO

[PopData] = TNC_ConvertTSDtoPopData('./',1);

%% Create continuous smoothed sdfs - also calculate low dimensional representation

currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 12;
% currParams.smthParams.decay    = 6;

% currParams.smthParams.decay    = 10;
currParams.filter.style = 0; % 0: Gaussian; 1: Causal Gaussian; 2: Boxcar

switch currParams.filter.style
    
    case 0
        [currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
        
    case 1
        [currParams.filter.kernel]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
        currParams.filter.kernel(1:currParams.smthParams.decay.*15) = 0;
                
    case 2
        currParams.filter.kernel = zeros(1,currParams.smthParams.decay.*3);
        currParams.filter.kernel(currParams.smthParams.decay:currParams.smthParams.decay.*2) = 1./currParams.smthParams.decay;

end


for ii=1:numel(PopData.session(1).unit)
    
    delta = zeros( 1 , round( size( ns4_data.Data , 2)./10 ) );
    if ii==1
        spk.density = zeros( numel(PopData.session(1).unit) , numel(delta) );
    end
        valids  = find(PopData.session(1).unit(ii).ts < numel(delta));
        tmp     = ceil( PopData.session(1).unit(ii).ts(valids) );
    delta( 1 , tmp ) = 1;
    
    spk.density(ii,:) = conv(delta,currParams.filter.kernel,'same');

end

[mappedA, mapping] = compute_mapping(spk.density', 'PCA', 3);

spk.pc_load = mappedA';
spk.pc_comp = mapping;

[v_pcl,i_pcl] = sort(spk.pc_comp.M(:,1));

%% TRAIN NN DECODER OF VELOCITY
clear train*
train_o=[]; train_i=[];

currParams.smthParams.rise     = 1;
currParams.smthParams.decay    = 10;
% currParams.smthParams.decay    = 6;

% currParams.smthParams.decay    = 10;
currParams.filter.style = 0; % 0: Gaussian; 1: Causal Gaussian; 2: Boxcar

switch currParams.filter.style
    
    case 0
        [currParams.filter.kernel2]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
        
    case 1
        [currParams.filter.kernel2]  = TNC_CreateGaussian(currParams.smthParams.decay.*15,currParams.smthParams.decay,currParams.smthParams.decay.*30,1);
        currParams.filter.kernel2(1:currParams.smthParams.decay.*15) = 0;
                
    case 2
        currParams.filter.kernel2 = zeros(1,currParams.smthParams.decay.*3);
        currParams.filter.kernel2(currParams.smthParams.decay:currParams.smthParams.decay.*2) = 1./currParams.smthParams.decay;

end

smth_vel = conv(reach.pos.v_win,currParams.filter.kernel2,'same');
smth_px = conv(reach.pos.xs,currParams.filter.kernel2,'same');
smth_py = conv(reach.pos.ys,currParams.filter.kernel2,'same');



% good_reaches = i_km_i(80:100);
% for mm=1:numel(good_reaches)
%     train_o = [train_o smth_vel(reach.start(good_reaches(mm)):reach.stop(good_reaches(mm)))];
%     train_i = [train_i spk.density(:,reach.start(good_reaches(mm)):reach.stop(good_reaches(mm)))];
% end
train_range = [84000:100000];
train_i = spk.density(:,train_range);
% train_i = spk.pc_load([1,3],train_range);
train_o = smth_vel(1,train_range);
% train_o = [smth_px(1,train_range) ; smth_py(1,train_range)]
figure(9); clf; plot(train_i(1,:).*1000); hold on; plot(train_o(1,:),'k','linewidth',2); plot(train_i(2,:).*1000); %plot(train_o(2,:),'r','linewidth',2);

%% TRAIN (SLOOOOW)

nn      = [33 12 6 1];
dIn     = [1];
dIntern = [1];
dOut    = [1];
net     = CreateNN(nn,dIn,dIntern,dOut);
netLM = train_LM( train_i , train_o , net, 10, 10);

decode_v = NNOut(train_i,netLM);

figure(21); clf;
plot(decode_v,'linewidth',2); hold on; plot(train_o,'k','linewidth',2);

%% FIND MOST SIMILAR REACHES AND MOST SIMIALR NEURAL RESPONSES AND COMPARE

for kk=1:numel(reach.start)
    
    spk.per_reach(kk,:) = spk.pc_load(1,reach.start(kk)-50:reach.start(kk)+250);
    mvmt.per_vel(kk,:)  = smth_vel(reach.start(kk)-50:reach.start(kk)+250);
end

% km_i = kmeans(spk.per_reach, 4);
% [v_km_i,i_km_i] = sort(km_i);

km_v = kmeans(mvmt.per_vel, 4);
[v_km_v,i_km_v] = sort(km_v);

% figure(10); clf; subplot(121); imagesc(spk.per_reach(i_km_i,:));
%  subplot(122); imagesc(corr(spk.per_reach(i_km_i,:)'));

figure(11); clf; subplot(221); imagesc(spk.per_reach(i_km_v,:));
 subplot(222); imagesc(corr(spk.per_reach(i_km_v,:)'));
 subplot(223); imagesc(mvmt.per_vel(i_km_v,:));
 subplot(224); imagesc(corr(mvmt.per_vel(i_km_v,:)'));

%% GRAB SOME REACH ALIGNED VIDEO TO MAKE HIGH QUALITY VIDEO EXAMPLES

frames.pre =25;
frames.post=35;
clear videos*;

valid_reach_starts = find( reach.start > min(ts_frames)+30000 & reach.start < max(ts_frames)-30000 );
disp(['Total reaches: ' num2str(numel(reach.start)) ' | Valid reaches: ' num2str(numel(valid_reach_starts))]);

[vals,inds] = sort(reach.stats.ampl(valid_reach_starts));

% find one example reach to look at
% example = [valid_reach_starts(inds(5)) valid_reach_starts(inds(2)) valid_reach_starts(inds(12)) valid_reach_starts(inds(10)) ]
% example = i_km_i([101,103,92,105])
example = i_km_i([60, 55, 58, 50]);
example = i_km_v([54, 55, 45, 61]);


for m=1:numel(example)
    example(m)
    video.key.ts        = double(reach.start(example(m))).*30;
    video.key.frame     = find(abs(ts_frames-video.key.ts)==min(abs(ts_frames-video.key.ts)),1);
    % find(ts_frames>video.key.ts,1)
    video.begin.frame   = video.key.frame - frames.pre;
    video.end.frame     = video.key.frame + frames.post;
    
    % LOAD THE VIDEO IN TO MEMORY
    videos{m}.vid = TNC_ReadSeqImages(seq_info, fid, [video.begin.frame:video.end.frame]);
    videos{m}.vid2= TNC_ReadSeqImages(seq_info2, fid2, [video.begin.frame:video.end.frame]);
    videos{m}.beg = video.begin.frame;
        begin1k   = reach.start(example(m))-(frames.pre*10);
        end1k     = reach.start(example(m))+(frames.post*10);
    videos{m}.r_v = smth_vel(begin1k:end1k);
    videos{m}.r_x = smth_px(begin1k:end1k);
    videos{m}.r_y = smth_py(begin1k:end1k);
    videos{m}.pcl = spk.pc_load(:,begin1k:end1k);    
    videos{m}.sdf = spk.density(:,begin1k:end1k);    
    videos{m}.d_v = NNOut(spk.density(:,begin1k:end1k),netLM);
%     videos{m}.d_v = NNOut(spk.pc_load([1,3],begin1k:end1k),netLM);
    videos{m}.d_v(1:20) = zeros(1,20);
    % LOAD THE VIDEO IN TO MEMORY
%     videos2{m}.vid = TNC_ReadSeqImages(seq_info2, fid2, [video.begin.frame:video.end.frame]);
    
end

%% Try to create a timelapse of images to see whether visualization works that way

tmp = medfilt2( videos{3}.vid{1}(:,:) ) - medfilt2( videos{3}.vid{1}(:,:) );

for j=2:numel([video.begin.frame:video.end.frame])-1
    tmp = tmp + medfilt2( videos{k}.vid{j}(:,:) ) - medfilt2( videos{3}.vid{j-1}(:,:) );
end

figure(100); imagesc(8268-tmp); colormap(gray);

%% PLAY VIDEO
cmap = TNC_CreateRBColormap(numel(i_pcl),'cpb');
figure(2); clf;
set(gcf,'Color','w');

v = VideoWriter('egReachPhys.avi');
v.FrameRate = 15;
v.Quality = 100;
open(v);

for j=2:numel([video.begin.frame:video.end.frame])-1
    for k=1:numel(videos)
        
        figure(2); subplot(5,4,k);
        movement    = medfilt2( videos{k}.vid2{j}(:,:) ) - medfilt2( videos{k}.vid2{j-1}(:,:) );
        structure   = medfilt2( videos{k}.vid2{j}(:,:) ) ;
        imagesc( movement + structure , [0 200] );
        title(videos{k}.beg+j); colormap(gray);
        axis off;
        figure(2); subplot(5,4,k+4);
        movement    = medfilt2( videos{k}.vid{j}(150:325,150:450) ) - medfilt2( videos{k}.vid{j-1}(150:325,150:450) );
        structure   = medfilt2( videos{k}.vid{j}(150:325,150:450) );
        imagesc( movement + structure , [0 200]);
        title(videos{k}.beg+j); colormap(gray);
        axis off;
        
        figure(2); subplot(5,4,k+8); hold off;
        %         plot(videos{k}.r_v(1:round(j.*10))); hold on;
        %         plot([j.*10 j.*10],[0 40],'r');
        %         axis([0 numel([video.begin.frame:video.end.frame]).*10 0 40]);
        plot([1:round(j.*10)]-400,(videos{k}.r_v(1:round(j.*10)).*10)-2500,'color',[0.7 0.7 0.7],'linewidth',2); hold on;
        plot(videos{k}.r_x(1:round(j.*10))-videos{k}.r_x(1),videos{k}.r_y(1:round(j.*10))-videos{k}.r_y(1),'linewidth',2); hold on;
        axis([-0.4e3 0.5e3 -2.5e3 0.5e3]);
        xlabel('x (a.u.) | Time (ms)');ylabel('y (a.u.)');
        
        figure(2); subplot(5,4,[k+12 k+16]); hold off;
        for jj=numel(i_pcl):-1:1
            patch( [0 1:round(j.*10) round(j.*10)] , [ 0 videos{k}.sdf( i_pcl(jj) , 1:round(j.*10 ) ) 0 ] + (jj.*0.025) , 'w' , 'EdgeColor' , cmap(i_pcl(jj),:) ); hold on;
            %             plot( videos{k}.sdf( i_pcl(jj) , 1:round(j.*10 ) ) + (jj.*0.025) , 'color' , cmap(i_pcl(jj),:) , 'linewidth' , 2 ); hold on;
        end
        %         plot(videos{k}.pcl(1,1:round(j.*10))-0.1, 'linewidth' , 2 ); hold on;
        %         plot(videos{k}.pcl(2,1:round(j.*10))-0.1, 'linewidth' , 2 ); hold on;
        plot((videos{k}.d_v(1,1:round(j.*10))./250)-0.5, 'k', 'linewidth' , 2 ); hold on;
        plot((videos{k}.r_v(1,1:round(j.*10))./250)-0.5, 'k', 'linewidth' , 1 , 'color' , [0.5 0.5 0.5]); hold on;
        %         plot(videos{k}.pcl(2,1:round(j.*10))); hold on;
        %         plot(videos{k}.pcl(3,1:round(j.*10))); hold on;
        %         plot([j.*10 j.*10],[-0.25 1],'r');
        axis([0 numel([video.begin.frame:video.end.frame]).*10 -0.5 1]);
        axis off;
        
    end
    drawnow;
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);


