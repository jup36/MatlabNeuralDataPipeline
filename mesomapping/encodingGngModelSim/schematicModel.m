%% Make toy networks to explore intuitions

% Selection network: separate population vectors for each movement direction / force
% Controller network: overlapping cosine tuned population vectors for each movement direction / force

% Compare 180 degree (~100%) separated actions (full dimensionality) and full force (1 - 10) range
% Compare 18 degree (10%) separated actions and 3g/6g (~25%) force range
clear D_*
pca_type = 'standard'
% pca_type = 'contrastive'


for iter = 1:10

[spk_map] = TNC_CreateRBColormap(8,'exag');
clear selector controller;
N = 120;
t_steps = 2500;
rew_t_gain = 500;
r_dir = [0 0 24 24];
r_tor = [3 6 3 6];
kernels = zeros(N*2,2e3);

t_lags_all = poissrnd(([1:N]*2));
tmp = sort(t_lags_all);
t_lags_reach = [tmp(1:2:end) tmp(2:2:end)] + 500;
t_lags_all = poissrnd(([1:N])*9);
tmp = sort(t_lags_all);
t_lags_pull = [tmp(1:2:end) tmp(2:2:end)] + 700;
figure(2); clf; subplot(121); plot(t_lags_reach); hold on; plot(t_lags_pull);

% Create random, variable temporal kernels
for ii=1:N
    [kernel_rise] = TNC_CreateGaussian(1e3,25+rand(1)*40,2e3,1);
    kernel_rise=kernel_rise/max(kernel_rise);
    [kernel] = TNC_CreateGaussian(1e3,50+rand(1)*200,2e3,1);
    kernel=kernel/max(kernel);
    kernel(1:1e3) = kernel_rise(1:1e3);
    kernels(ii,:) = kernel;

    [kernel_rise] = TNC_CreateGaussian(1e3,rand(1)*(t_lags_pull(ii)-600)/4,2e3,1);
    kernel_rise=kernel_rise/max(kernel_rise);
    [kernel] = TNC_CreateGaussian(1e3,50+rand(1)*300,2e3,1);
    kernel=kernel/max(kernel);
    kernel(1:1e3) = kernel_rise(1:1e3);
    kernels(ii+N,:) = kernel;
end
figure(2); subplot(122); imagesc(kernels);


% For meaningful control encoding scheme set gains
gains = rand(1,N).*0.8;
phases = -pi:pi/N:pi;
rew_gains = poissrnd(3,1,N).*(gains);

for jj=1:numel(r_dir)

    for ii=1:N

    % build separate population vectors
        delta = zeros(1,t_steps);
        type = find(r_dir(jj)==unique(r_dir));
        if type>2
            type=2;
        end
        if ii<(type*60) & ii>((type-1)*60)
            trial_amp_reach = randperm(10,1);
        else
            trial_amp_reach = rand(1)*2;
        end
        delta(t_lags_reach(ii)) = trial_amp_reach;

        % for all trials add a probabilistic response to reward
        delta(round(1750+rand(1)*rew_t_gain)) = rew_gains(ii);

        selector.psth(ii,:,jj) = conv( delta , kernels(ii,:) , 'same' );

        delta = zeros(1,t_steps);
        type = find(r_tor(jj)==unique(r_tor));
        if type>2
            type=2;
        end
        if ii<(type*60) & ii>((type-1)*60)
            trial_amp_pull = randperm(8,1);
        else
            trial_amp_pull = rand(1)*2;
        end
        delta(t_lags_pull(ii)) = trial_amp_pull;

        % for all trials add a probabilistic response to reward
        delta(round(1750+rand(1)*rew_t_gain)) = rew_gains(ii);

        selector.psth(ii+N,:,jj) = conv( delta , kernels(ii+N,:) , 'same' );


    % build cosine tuned and monotonic force tuned psths
        delta = zeros(1,t_steps);
        angle_in_radians = r_dir(jj)/180*pi;
        trial_amp_reach = randperm(10,1)*cos( (phases(ii)+angle_in_radians)/1.75 );
        if trial_amp_reach<0
            trial_amp_reach=0;
        end
        delta(t_lags_reach(ii)) = trial_amp_reach;

        % for all trials add a probabilistic response to reward
        delta(round(1750+rand(1)*rew_t_gain)) = rew_gains(ii);

        controller.psth(ii,:,jj) = conv( delta , kernels(ii,:) , 'same' );

        delta = zeros(1,t_steps);
        trial_gn_pull = cos( (phases(ii)+(angle_in_radians/2))/1.75 );
        trial_amp_pull = gains(ii)*r_tor(jj)*trial_gn_pull;
        delta(t_lags_pull(ii)) = trial_amp_pull;

        % for all trials add a probabilistic response to reward
        delta(round(1750+rand(1)*rew_t_gain)) = rew_gains(ii);

        controller.psth(ii+N,:,jj) = conv( delta , kernels(ii+N,:) , 'same' );

    end
end

% Plotting output of networks in full D and PC projections

figure(3); clf;

% fake_depth_inds = randperm(size(selector.psth,1));
% fake_depth_inds = 1:size(selector.psth,1);
[~,tmp] = max(mean(selector.psth(:,1:1650,:),3),[],2);
[~,fake_depth_inds] = sort(tmp,'ascend');

dirs = unique(r_dir);
for kk=dirs
    subplot(2,numel(dirs),find(kk==dirs))
    imagesc(mean(selector.psth(fake_depth_inds,:,r_dir==kk),3),[0 10]);
    colormap(spk_map); box off; ylabel('Unit index'); xlabel('timesteps');

    subplot(2,numel(dirs),find(kk==dirs)+numel(dirs))
    imagesc(mean(controller.psth(fake_depth_inds,:,r_dir==kk),3),[0 10]);
    colormap(spk_map); box off; ylabel('Unit index'); xlabel('timesteps');
end

figure(10); clf;
plot(-500:1999,mean(mean(controller.psth(fake_depth_inds,:,:),3),1),'LineWidth',2); hold on;
plot(-500:1999,mean(mean(selector.psth(fake_depth_inds,:,:),3),1),'LineWidth',2); hold on;
legend({'Controller','Selector'});  xlabel('timesteps'); ylabel('Mean activity'); box off;

figure(4); clf;
tors = unique(r_tor);
for kk=tors
    subplot(2,numel(tors),find(kk==tors))
    imagesc(mean(selector.psth(fake_depth_inds,:,r_tor==kk),3),[0 10]);

    subplot(2,numel(tors),find(kk==tors)+numel(tors))
    imagesc(mean(controller.psth(fake_depth_inds,:,r_tor==kk),3),[0 10]);
end

trialgroup.spec(iter).psth = controller.psth;
trialgroup.selt(iter).psth = selector.psth;

end



%% Plotting routine

switch pca_type
    case 'standard'
        % Consider projection onto PCs for experiment with big difference in
        % direction
        [mA_sel, m_sel] = compute_mapping(reshape(selector.psth(:,:,[1 3 7 9]),[],4*t_steps)','PCA',3);
        [mA_con, m_con] = compute_mapping(reshape(controller.psth(:,:,[1 3 7 9]),[],4*t_steps)','PCA',3);
        
        % Consider projection onto PCs for experiment with small difference in
        % direction & torque
        [mA_sel_jp, m_sel_jp] = compute_mapping(reshape(selector.psth(:,:,[1 2 4 5]),[],4*t_steps)','PCA',3);
        [mA_con_jp, m_con_jp] = compute_mapping(reshape(controller.psth(:,:,[1 2 4 5]),[],4*t_steps)','PCA',3);

    case 'contrastive'
        % Consider projection onto PCs for experiment with big difference in
        % direction
        [mA_sel, m_sel] = compute_mapping([reshape(selector.psth(:,:,3)-selector.psth(:,:,1),[],t_steps) reshape(selector.psth(:,:,9)-selector.psth(:,:,7),[],t_steps)]','PCA',3);
        [mA_con, m_con] = compute_mapping([reshape(controller.psth(:,:,3)-controller.psth(:,:,1),[],t_steps) reshape(controller.psth(:,:,9)-controller.psth(:,:,7),[],t_steps)]','PCA',3);
        
        % Consider projection onto PCs for experiment with small difference in
        % direction & torque
        [mA_sel_jp, m_sel_jp] = compute_mapping([reshape(selector.psth(:,:,2)-selector.psth(:,:,1),[],t_steps) reshape(selector.psth(:,:,5)-selector.psth(:,:,4),[],t_steps)]','PCA',3);
        [mA_con_jp, m_con_jp] = compute_mapping([reshape(controller.psth(:,:,2)-controller.psth(:,:,1),[],t_steps) reshape(controller.psth(:,:,5)-controller.psth(:,:,4),[],t_steps)]','PCA',3);

end

catmap      = [ 58 84 162 255; 233 41 42 255 ; 91 116 178 255; 229 136 140 255]/255;

ctx_str_cnxn_type = 'random'

% common_view = [-60 20];
% common_view = [-28 -22]   ;
% common_view = [-120 70];
common_view = [-110 20]
common_scaling_sel = [-40 40 -25 35 -25 30];
common_scaling_con = [-60 20 -35 5 -25 25];

figure(5); clf;
subplot(221);
% plot PC projection for 4 trial conditions in selector big difference
% dir Left tor Low
plot3(selector.psth(:,:,1)'*m_sel.M(:,1),selector.psth(:,:,1)'*m_sel.M(:,2),selector.psth(:,:,1)'*m_sel.M(:,3),'color',catmap(3,:),'LineWidth',2); hold on;
plot3(selector.psth(:,:,3)'*m_sel.M(:,1),selector.psth(:,:,3)'*m_sel.M(:,2),selector.psth(:,:,3)'*m_sel.M(:,3),'color',catmap(1,:),'LineWidth',2);
plot3(selector.psth(:,:,7)'*m_sel.M(:,1),selector.psth(:,:,7)'*m_sel.M(:,2),selector.psth(:,:,7)'*m_sel.M(:,3),'color',catmap(4,:),'LineWidth',2); hold on;
plot3(selector.psth(:,:,9)'*m_sel.M(:,1),selector.psth(:,:,9)'*m_sel.M(:,2),selector.psth(:,:,9)'*m_sel.M(:,3),'color',catmap(2,:),'LineWidth',2);
view(common_view);
% axis(common_scaling_sel);
grid on; xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
tmp_big_diff = procrustes([selector.psth(:,:,3)'*m_sel.M(:,1),selector.psth(:,:,3)'*m_sel.M(:,2),selector.psth(:,:,3)'*m_sel.M(:,3)]', [selector.psth(:,:,9)'*m_sel.M(:,1),selector.psth(:,:,9)'*m_sel.M(:,2),selector.psth(:,:,9)'*m_sel.M(:,3)]','Scaling',0,'reflection',0)
title(['Selector | Large $\Delta$ | ' num2str(tmp_big_diff)],'Interpreter','Latex');

subplot(222);
% plot PC projection for 4 trial conditions in controller jp difference
plot3(controller.psth(:,:,1)'*m_con.M(:,1),controller.psth(:,:,1)'*m_con.M(:,2),controller.psth(:,:,1)'*m_con.M(:,3),'color',catmap(3,:),'LineWidth',2); hold on;
plot3(controller.psth(:,:,3)'*m_con.M(:,1),controller.psth(:,:,3)'*m_con.M(:,2),controller.psth(:,:,3)'*m_con.M(:,3),'color',catmap(1,:),'LineWidth',2);
plot3(controller.psth(:,:,7)'*m_con.M(:,1),controller.psth(:,:,7)'*m_con.M(:,2),controller.psth(:,:,7)'*m_con.M(:,3),'color',catmap(4,:),'LineWidth',2); hold on;
plot3(controller.psth(:,:,9)'*m_con.M(:,1),controller.psth(:,:,9)'*m_con.M(:,2),controller.psth(:,:,9)'*m_con.M(:,3),'color',catmap(2,:),'LineWidth',2);
view(common_view);
% axis(common_scaling_con);
grid on; xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
tmp_small_diff = procrustes([controller.psth(:,:,3)'*m_con.M(:,1),controller.psth(:,:,3)'*m_con.M(:,2),controller.psth(:,:,3)'*m_con.M(:,3)]', [controller.psth(:,:,9)'*m_con.M(:,1),controller.psth(:,:,9)'*m_con.M(:,2),controller.psth(:,:,9)'*m_con.M(:,3)]','Scaling',0,'reflection',0)
title(['Controller | Large $\Delta$ | ' num2str(tmp_small_diff)],'Interpreter','Latex');


subplot(223);
% plot PC projection for 4 trial conditions in selector
plot3(selector.psth(:,:,1)'*m_sel_jp.M(:,1),selector.psth(:,:,1)'*m_sel_jp.M(:,2),selector.psth(:,:,1)'*m_sel_jp.M(:,3),'color',catmap(3,:),'LineWidth',2); hold on;
plot3(selector.psth(:,:,2)'*m_sel_jp.M(:,1),selector.psth(:,:,2)'*m_sel_jp.M(:,2),selector.psth(:,:,2)'*m_sel_jp.M(:,3),'color',catmap(1,:),'LineWidth',2);
plot3(selector.psth(:,:,4)'*m_sel_jp.M(:,1),selector.psth(:,:,4)'*m_sel_jp.M(:,2),selector.psth(:,:,4)'*m_sel_jp.M(:,3),'color',catmap(4,:),'LineWidth',2); hold on;
plot3(selector.psth(:,:,5)'*m_sel_jp.M(:,1),selector.psth(:,:,5)'*m_sel_jp.M(:,2),selector.psth(:,:,5)'*m_sel_jp.M(:,3),'color',catmap(2,:),'LineWidth',2);
view(common_view);
% axis(common_scaling_sel);
grid on; xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
D_big_diff(iter) = procrustes([selector.psth(:,:,2)'*m_sel_jp.M(:,1),selector.psth(:,:,2)'*m_sel_jp.M(:,2),selector.psth(:,:,2)'*m_sel_jp.M(:,3)]', [selector.psth(:,:,5)'*m_sel_jp.M(:,1),selector.psth(:,:,5)'*m_sel_jp.M(:,2),selector.psth(:,:,5)'*m_sel_jp.M(:,3)]','Scaling',0,'reflection',0)
title(['Selector | Small $\Delta$ | ' num2str(D_big_diff(iter))],'Interpreter','Latex');


subplot(224);
% plot PC projection for 4 trial conditions in controller big difference
plot3(controller.psth(:,:,1)'*m_con_jp.M(:,1),controller.psth(:,:,1)'*m_con_jp.M(:,2),controller.psth(:,:,1)'*m_con_jp.M(:,3),'color',catmap(3,:),'LineWidth',2); hold on;
plot3(controller.psth(:,:,2)'*m_con_jp.M(:,1),controller.psth(:,:,2)'*m_con_jp.M(:,2),controller.psth(:,:,2)'*m_con_jp.M(:,3),'color',catmap(1,:),'LineWidth',2);
plot3(controller.psth(:,:,4)'*m_con_jp.M(:,1),controller.psth(:,:,4)'*m_con_jp.M(:,2),controller.psth(:,:,4)'*m_con_jp.M(:,3),'color',catmap(4,:),'LineWidth',2); hold on;
plot3(controller.psth(:,:,5)'*m_con_jp.M(:,1),controller.psth(:,:,5)'*m_con_jp.M(:,2),controller.psth(:,:,5)'*m_con_jp.M(:,3),'color',catmap(2,:),'LineWidth',2);
view(common_view);
% axis(common_scaling_con);
grid on; xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
D_small_diff(iter) = procrustes([controller.psth(:,:,2)'*m_con_jp.M(:,1),controller.psth(:,:,2)'*m_con_jp.M(:,2),controller.psth(:,:,2)'*m_con_jp.M(:,3)]', [controller.psth(:,:,5)'*m_con_jp.M(:,1),controller.psth(:,:,5)'*m_con_jp.M(:,2),controller.psth(:,:,5)'*m_con_jp.M(:,3)]','Scaling',0,'reflection',0)
title(['Controller | Small $\Delta$ | ' num2str(D_small_diff(iter))],'Interpreter','Latex');

print([pca_type '/Model-PCs-' num2str(iter)],'-dsvg','-painters')



save([pca_type '/RotationVals'],'D_*');

figure(100); clf;
boxplot([(D_big_diff)',(D_small_diff)'],{'Selection','Specification'});
ylabel('Intertarget Rotation (norm.)'); box off;
print('Model-PCs-SummaryDistance','-dsvg','-painters');
title(pca_type); box off;

%% Plot the data in the same fashion for comparing mean PSTH and sorted PSTHs

% get cortex psths and then sort by latency
ctx_psths = mtv.all_psth(mtv.all_labels==0,:);
[~,lags] = max(ctx_psths,[],2);
[~,cis] = sort(lags,'ascend');

% get str psths and then sort by latency
str_psths = mtv.all_psth(mtv.all_labels==1,:);
[~,lags_str] = max(str_psths,[],2);
[~,sis] = sort(lags_str,'ascend');

% concatenate response
figure(125); clf;
imagesc([ ctx_psths(cis,:) ;  str_psths(sis,:) ],[0 4]); colormap(spk_map);
% imagesc(str_psths(sis,:),[0 5]); colormap(spk_map);
ylabel('Sorted units'); xlabel('Time from pull (ms)');
axis([1250 3750 0 size([ ctx_psths(cis,:) ;  str_psths(sis,:) ],1)]);
% axis([1250 3750 0 size(str_psths,1)]);
box off;

figure(126); clf;
plot(mean([ ctx_psths(cis,:) ;  str_psths(sis,:) ],1),'linewidth',2);
axis([1250 3750 -0.15 0.75]); 
ylabel('Mean activity (z)'); xlabel('Time from pull (ms)'); box off;

%% Conceputal visualization of tuning either on cosine/control model or high D projection

pref_dir = [150:-50:-150]*pi/360;
angles = -pi/1.5:0.01:pi/1.5;
clear rad;

figure(300); clf;
figure(301); clf;

even_dir = [150:-50:-150]*pi/180;

rand_ths = even_dir(randperm(numel(even_dir)));
rand_phis = even_dir(randperm(numel(even_dir)));

colorscheme = [72 140 203;
    108 204 222;
    109 194 128;
    155 203 60;
    248 236 32;
    247 147 30;
    238 50 36] / 256;


for gg=1:numel(pref_dir)
   
    rad(gg,:) = cos((angles-pref_dir(gg))).^3;
    rad(gg,rad(gg,:)<=0) = 0;
%     polarplot(angles,rad(gg,:),'color',colorscheme(gg,:),'linewidth',2); hold on;
    
figure(300);
    [x,y,z] = sph2cart(angles,zeros(1,numel(angles)),rad(gg,:));    
    plot3(x,y,z,'color',colorscheme(gg,:),'linewidth',2); hold on;    
    box off; grid on; axis([-1 1 -1 1 -1 1]);
    
figure(301);
    [x,y,z] = sph2cart([rand_ths(gg) rand_ths(gg)],[rand_phis(gg) rand_phis(gg)],[0 1]);    
    plot3(x,y,z,'color',colorscheme(gg,:),'linewidth',2); hold on;    
    box off; grid on; axis([-1 1 -1 1 -1 1]);
    
end

%% Take the encoding model and turn it into a network and ask whether decoding analysis recovers the W matrix
[spk_map] = TNC_CreateRBColormap(8,'bb-sym');

str_as_ctx_driven = 0;
str_ctx_noise = 0;
scaler = 1.45;
regen = 1;

if regen
    
    gains_tmp = rand(50,1)-0.15;
    gains_sort = sort(gains_tmp,'descend');

    str_inds = randperm(50);
    gains_tmp2 = rand(50,1)-0.1;
    gains_sort2 = sort(gains_tmp2,'descend');
    gains_sort_str = gains_sort2(str_inds);
    
    weights = sort(rand(1,50)-0.1,'descend');
%     gains_sort_str = sort(gains_tmp,'ascend');

    move_noise = randn(1,100);
    move_noiseS = randn(1,100);
    
    for kk=1:50        
        % io curves with distinct gains I guess is a decent model?
        ctx(kk,:) = scaler.*randn(1,100)+(gains_sort(kk).*[1:100]./15);
        str_ind(kk,:) = scaler.*randn(1,100)+(gains_sort_str(kk).*[1:100]./15);
    end

%     W_ctx = gains_sort';
%     W_str = gains_sort_str';

%     weights = poissrnd(1,1,50);
%     weights = weights ./ max(weights);
%     W_ctx = sort(weights,'descend');
%     W_str = W_ctx(str_inds);
    
    W_ctx = sort(rand(1,50)-0.1,'descend');
    W_str = weights(str_inds);
%     W_str = sort(rand(1,50)-0.1,'ascend');

%     W_ctx = rand(1,50)-0.1;
%     W_str = rand(1,50)-0.1;

end

ctx_str_cnxn_type = 'distrib'

% create a biased local connectivity matrix
switch ctx_str_cnxn_type

    case 'cent_surr'
        W_ctx_str = diag(ones(1,50));
        cnxn_kern_e = TNC_CreateGaussian(25,2,50,1);
        cnxn_kern_i = TNC_CreateGaussian(25,16,50,1);
        cnxn_kern = cnxn_kern_e./max(cnxn_kern_e) - 0.5*cnxn_kern_i./max(cnxn_kern_i);
        for kk=1:50
            W_ctx_str(kk,:) = conv(W_ctx_str(kk,:),cnxn_kern,'same');
        end

    case 'cent_surr_rand'
        W_ctx_str = diag(ones(1,50));
        cnxn_kern_e = TNC_CreateGaussian(25,2,50,1);
        cnxn_kern_i = TNC_CreateGaussian(25,10,50,1);
        cnxn_kern = cnxn_kern_e./max(cnxn_kern_e) - 0.5*cnxn_kern_i./max(cnxn_kern_i);
        for kk=1:50
            W_ctx_str(kk,:) = conv(W_ctx_str(kk,:),cnxn_kern,'same');
        end
        W_ctx_str = W_ctx_str - (rand(50,50)-0.5);
        
    case 'topo_rand'
        W_ctx_str = diag(ones(1,50));
        cnxn_kern = TNC_CreateGaussian(25,2,50,1);
        cnxn_kern = cnxn_kern./max(cnxn_kern);
        for kk=1:50
            W_ctx_str(kk,:) = conv(W_ctx_str(kk,:),cnxn_kern,'same');
        end
        W_ctx_str = W_ctx_str - (rand(50,50)-0.5);

    case 'distrib'
        W_ctx_str = diag(ones(1,50));
        cnxn_kern_e(1,:)  = TNC_CreateGaussian(25,2,50,1);
        cnxn_kern_e(2,:)  = TNC_CreateGaussian(37,3,50,1)/2;
        cnxn_kern_e(3,:)  = TNC_CreateGaussian(49,3,50,1)/3;
        cnxn_kern_e(4,:)  = TNC_CreateGaussian(13,4,50,1)/2;
        cnxn_kern_e(5,:)  = TNC_CreateGaussian(1,4,50,1)/3;
%         cnxn_kern = sum(cnxn_kern_e,1)./max(sum(cnxn_kern_e,1));

        cnxn_kern_i(1,:)  = TNC_CreateGaussian(25,10,50,1);
        cnxn_kern_i(2,:)  = TNC_CreateGaussian(37,10,50,1)/2;
        cnxn_kern_i(3,:)  = TNC_CreateGaussian(49,10,50,1)/3;
        cnxn_kern_i(4,:)  = TNC_CreateGaussian(13,10,50,1)/2;
        cnxn_kern_i(5,:)  = TNC_CreateGaussian(1,10,50,1)/3;
        cnxn_kern = sum(cnxn_kern_e,1)./max(sum(cnxn_kern_e,1)) - 0.25*sum(cnxn_kern_i,1)./max(sum(cnxn_kern_i,1));
        for kk=1:50
            W_ctx_str(kk,:) = conv(W_ctx_str(kk,:),cnxn_kern,'same');
        end
        W_ctx_str = W_ctx_str - 0.2 - (rand(50,50)-0.5);
        
    case 'random'
        W_ctx_str = rand(50,50)-0.5;

    case 'topographic'
        W_ctx_str = diag(ones(1,50));
        cnxn_kern = TNC_CreateGaussian(25,10,50,1);
        cnxn_kern = cnxn_kern./max(cnxn_kern);
        for kk=1:50
            W_ctx_str(kk,:) = conv(W_ctx_str(kk,:),cnxn_kern,'same');
        end
        W_ctx_str = W_ctx_str;
end

if str_as_ctx_driven == 1
    str = (W_ctx_str*ctx);
elseif str_as_ctx_driven == -1
    str = str_ind;
else
    str = 0.9.*str_ind + 0.1.*(W_ctx_str*ctx);
end

if str_ctx_noise

    switchy = randperm(100);
    str(:,switchy(1:25)) = str(:,switchy(1:25))./scaler; 
    ctx(:,switchy(1:25)) = ctx(:,switchy(1:25)).*scaler;
    switchy = randperm(100);
    str(:,switchy(51:100)) = str(:,switchy(51:100)).*scaler; 
    ctx(:,switchy(51:100)) = ctx(:,switchy(51:100))./scaler; 
    
end

Cmop = triu( corr(ctx') , 1 );
Cstr = triu( corr(str') , 1 );

figure(99); 
subplot(151); imagesc(ctx); title('ctx');
subplot(152); imagesc(str); title('str');
subplot(153); imagesc(W_ctx_str,[-1 1]); colormap(spk_map); title('W_{c,s}');
subplot(154); imagesc(corr(ctx'),[-1 1]); colormap(spk_map); title(['corr(ctx) : ' num2str( mean(abs(Cmop(Cmop~=0))) ) ]);
subplot(155); imagesc(corr(str'),[-1 1]); colormap(spk_map); title(['corr(str) : ' num2str( mean(abs(Cstr(Cstr~=0))) ) ]);


%% Examine RRR for the model

% compute PCs of the W_ctx_str matrix
[COEFF, SCORE] = pca(W_ctx_str);
str = 0.8.*str_ind + 0.2.*(W_ctx_str*ctx);
% str = (W_ctx_str*ctx);

figure(100); clf;
subplot(131); imagesc(COEFF); 
subplot(132); imagesc(SCORE); 
subplot(133); imagesc(SCORE*COEFF'); 
colormap(spk_map);

figure(101); clf;
num_pcs = [1 2 5 10 20 49 ];
r_sq = [];
for kk=1:6
    tmp_COEFF = COEFF;
    tmp_COEFF(:,num_pcs(kk)+1:end) = 0;
    subplot(3,6,kk);
    imagesc(tmp_COEFF);
    
    subplot(3,6,kk+6);
    imagesc(SCORE*tmp_COEFF');
    colormap(spk_map);

    subplot(3,6,kk+12);
    str_rrr = 0.8.*str_ind + 0.2.*((SCORE*tmp_COEFF')*ctx);
    % str_rrr = (SCORE*tmp_COEFF')*ctx;
    imagesc(str_rrr);
    title(num2str( (corr2(str,str_rrr)).^2) );

    r_sq(kk) = (corr2(str,str_rrr)).^2;
end

r_sq = [];
for kk=1:49
    tmp_COEFF = COEFF;
    tmp_COEFF(:,kk+1:end) = 0;
    str_rrr = (SCORE*tmp_COEFF')*ctx;
    r_sq(kk) = (corr2(str,str_rrr)).^2;
end

figure(102); clf;
scatter(1:49,r_sq,50,'k','filled');

% Examine how much variance can explain of observed str matrix based upon 1:N PCs from W_ctx_str

%% Different conceptual models of how movement is generated

% ctx_decode_inds = randperm(50,25);
% str_decode_inds = randperm(50,25);
regenM=0;
Mscaler = 15;
if regenM
    noise = Mscaler.*randn(100,1);
end

% biasing model (sums downstream)
model_id = 'equal sum'
movement = (str'*W_str') + (ctx'*W_ctx') + noise;
% 
% model_id = 'str bias';
% movement = 1.33*(str'*W_str') + 0.67*(ctx'*W_ctx') + noise;

% model_id = 'str only'
% movement = 2*(str'*W_str') +noise;
 
% model_id = 'ctx bias'
% movement = 0.67*(str'*W_str') + 1.33*(ctx'*W_ctx') + noise;

% model_id = 'ctx only'
% movement = 2*(ctx'*W_ctx') + noise;

% reconstruct
committee = 2; batch=200;
clear *_tmp
[cat_map] = TNC_CreateRBColormap(8,'cat2');

% movement = movement + 0.2.*movement.*randn(100,1);


switch committee

    case 1
        for qq=1:batch
            d_ts = randperm(90,30)+5;
            W_hat_str_tmp(qq,:) = pinv(str(:,d_ts)')*movement(d_ts);
            W_hat_ctx_tmp(qq,:) = pinv(ctx(:,d_ts)')*movement(d_ts);
        end
        W_hat_str = mean(W_hat_str_tmp,1)';
        W_hat_ctx = mean(W_hat_ctx_tmp,1)';
        W_hat_all = [W_hat_str ; W_hat_ctx];

    case 2
        for qq=1:batch
            d_ts = randperm(90,30)+5;
            W_hat_all_tmp(qq,:) = pinv([str(:,d_ts) ; ctx(:,d_ts)]')*movement(d_ts);
        end
        W_hat_all = mean(W_hat_all_tmp,1)';
        W_hat_str = W_hat_all(1:50);
        W_hat_ctx = W_hat_all(51:100);
        
    case 0
        W_hat_all = pinv([str ; ctx]')*movement;
        W_hat_str = W_hat_all(1:50);
        W_hat_ctx = W_hat_all(51:100);

    case -1
        W_hat_str = pinv(str')*movement;
        W_hat_ctx = pinv(ctx')*movement;
        
end

% movement_hat = [str ; ctx]'*W_hat_all;
movement_hat = (str'*W_hat_str) + (ctx'*W_hat_ctx);

figure('Name',[model_id '-' ctx_str_cnxn_type]); clf;

subplot(2,4,[1 2 5 6]);
plot([-50 250],[-50 250],'k-'); hold on;
if committee>0
    scatter(movement,movement_hat,75,[zeros(1,5) ones(1,90) zeros(1,5)],'filled'); box off; colormap(cat_map([1 4 5 6],:)); hold on;
    scatter(movement,movement_hat-(str'*W_hat_str),40,[zeros(1,5) ones(1,90) zeros(1,5)]+2,'filled');
else
    scatter(movement,movement_hat,60,'filled'); box off; colormap(cat_map([1 4],:))
end
ylabel('Predicted')
xlabel('Observed')
title('Full model decode')
axis([-50 250 -50 250]);

subplot(2,4,3);
plot([-50 250],[-50 250],'k--'); hold on;
scatter(movement,ctx'*W_hat_ctx,'filled');
m1_rmse = sqrt( mean( (movement - ctx'*W_hat_ctx).^2 ) );
PT = polyfit(movement,ctx'*W_hat_ctx,1);
title(['M1 tuning: ' num2str(PT(1),'%.2f')]);
ylabel('Predicted')
xlabel('Observed')
axis([-50 250 -50 250]);



subplot(2,4,4);
plot([-50 250],[-50 250],'k--'); hold on;
scatter(movement,str'*W_hat_str,'filled');
str_rmse = sqrt( mean( (movement - str'*W_hat_str).^2 ) );
PT = polyfit(movement,str'*W_hat_str,1);
title(['STR tuning: ' num2str(PT(1),'%.2f')]);
ylabel('Predicted')
xlabel('Observed')
axis([-50 250 -50 250]);



subplot(2,4,7);
plot([-0.5 2],[-0.5 2],'k--'); hold on;
scatter(W_ctx,W_hat_ctx);
axis([-0.5 1 -0.5 1]);
title('Inferred W_{ctx}')
ylabel('Predicted')
xlabel('Observed')

subplot(2,4,8);
plot([-0.5 2],[-0.5 2],'k--'); hold on;
scatter(W_str,W_hat_str);
axis([-0.5 1 -0.5 1]);
title('Inferred W_{str}')
ylabel('Predicted')
xlabel('Observed')

% set(gcf,'OuterPosition',[871 1539-(10*get(gcf,'Number')) 838 499]);


%% Look at data and see what STR vs CTX heatmaps of trialwise response look like
sess = 4;
[exag_map] = TNC_CreateRBColormap(8,'exag');
[sign_map] = TNC_CreateRBColormap(8,'bb-sym');

[fs,fis] = sort(mtv.stamps(sess).max_f,'ascend');


for kk=1:size(mtv.stamps(sess).pullAct,1)
    P = polyfit(mtv.stamps(sess).max_f,mtv.stamps(sess).pullAct(kk,:),1);
    mtv.stamps(sess).pullAct_slopes(kk) = P(1);
end

m1.slopes = mtv.stamps(sess).pullAct_slopes(mtv.stamps(sess).spkTimesID==0);
m1.pullAct = mtv.stamps(sess).pullAct(mtv.stamps(sess).spkTimesID==0,:);
[m1.ss,m1.sis] = sort(m1.slopes,'descend');

str.slopes = mtv.stamps(sess).pullAct_slopes(mtv.stamps(sess).spkTimesID==1);
str.pullAct = mtv.stamps(sess).pullAct(mtv.stamps(sess).spkTimesID==1,:);
[str.ss,str.sis] = sort(str.slopes,'descend');


Cmop = triu( corr(m1.pullAct(m1.sis,fis)') , 1 );
Cstr = triu( corr(str.pullAct(str.sis,fis)') , 1 );


figure(100); clf;
subplot(141);
imagesc(m1.pullAct(m1.sis,fis),[-8000 10000]);
subplot(142);
imagesc(str.pullAct(str.sis,fis),[-8000 10000]);
subplot(143);
imagesc(corr(m1.pullAct(m1.sis,fis)')),[-1 1]; 
title(['corr(ctx) : ' num2str( mean(abs(Cmop(Cmop~=0))) ) ]);
subplot(144);
imagesc(corr(str.pullAct(str.sis,fis)'),[-1 1]); 
title(['corr(str) : ' num2str( mean(abs(Cstr(Cstr~=0))) ) ]);
colormap(sign_map);

% figure(101); clf;
% plot(fs,sum(m1.pullAct,1)); hold on; plot(fs,sum(str.pullAct,1));

%% Use this to make a figure 1 schematic

%% Load pkl data from Catia into matlab


fid=py.open(filename,'rb');
data=py.pickle.load(fid);

% or 

pickle = py.importlib.import_module('pickle');
fh = py.open('data.pkl', 'rb')
P = pickle.load(fh);    % pickle file loaded to Python variable
fh.close();