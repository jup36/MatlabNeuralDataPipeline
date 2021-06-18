
np11_pullStarts = load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/WR39_100219/binSpkCountSTRWR39_100219.mat'),'pullStarts'); 
np11_pullStarts = np11_pullStarts.pullStarts; 

np12_pullStarts = load(fullfile('/Volumes/Beefcake/Junchol_Data/JS2p0/WR40_081919/Matfiles/binSpkCountSTRWR40_081919.mat'),'pullStarts'); 
np12_pullStarts = np12_pullStarts.pullStarts; 

np2_pullStarts = load(fullfile('/Volumes/Beefcake/Junchol_Data/np20_test/WR39_100819_g0/WR39_100819_g0_imec0/binSpkCountWR39_100819.mat'),'pullStarts'); 
np2_pullStarts = np2_pullStarts.pullStarts; 

clearvars -except np11_pullStarts np12_pullStarts np2_pullStarts

%% Compare single chank and 4 shank data 
spk_kernel = TNC_CreateGaussian(500,40,1000,1);
[lat_map] = TNC_CreateRBColormap(4000,'cpb');
[cat_map] = TNC_CreateRBColormap(4,'mapb');
[phys_map] = TNC_CreateRBColormap(20,'mbr');
for k=1:numel(np2_pullStarts.geometry)    
    np2_pullStarts.xpos(k) = 250*floor(np2_pullStarts.geometry{k}(1)/250);
    np2_pullStarts.ypos(k) = round(np2_pullStarts.geometry{k}(2));    
end
for k=1:numel(np11_pullStarts.geometry)
    np11_pullStarts.xpos(k) = 250*floor(np11_pullStarts.geometry{k}(1)/250);
    np11_pullStarts.ypos(k) = round(np11_pullStarts.geometry{k}(2));    
end
for k=1:numel(np12_pullStarts.geometry)
    np12_pullStarts.xpos(k) = 250*floor(np12_pullStarts.geometry{k}(1)/250);
    np12_pullStarts.ypos(k) = round(np12_pullStarts.geometry{k}(2));    
end
% only going to consider dorsal striatum units
np11_pullStarts.valid_units = find(np11_pullStarts.ypos>2.2e3 & np11_pullStarts.ypos<3.1e3);    
np12_pullStarts.valid_units = find(np12_pullStarts.ypos>2.2e3 & np12_pullStarts.ypos<3.1e3);    
np2_pullStarts.valid_units = find(np2_pullStarts.ypos>2.2e3 & np2_pullStarts.ypos<3.15e3);    
% create concatenated sdf for computing PCA

% 1-shank 1st session
pts = 1001:5000; %-- UNIQUE TO THIS DATASET; CHANGE BACK LATER (-3000 TO 2000 ms)
np11_struct.cat = []; 
for jj=1:size(np11_pullStarts.unitTimeTrial,3)
    spk_mat = zeros(numel(np11_pullStarts.valid_units),numel(pts));
    
    for mm=1:numel(np11_pullStarts.valid_units)
        if jj==1
            np11_struct.lag(mm) = find(np11_pullStarts.trAVG{np11_pullStarts.valid_units(mm)}(1,pts)==max(np11_pullStarts.trAVG{np11_pullStarts.valid_units(mm)}(1,pts)),1);
            np11_struct.mag(mm) = max(np11_pullStarts.trAVG{np11_pullStarts.valid_units(mm)}(1,pts))-mean(np11_pullStarts.trAVG{np11_pullStarts.valid_units(mm)}(1,pts));
        end
        spk_mat(mm,:) = conv( np11_pullStarts.unitTimeTrial(np11_pullStarts.valid_units(mm),pts,jj) ,spk_kernel,'same');   
    end
    
    np11_struct.wins(jj).win_dat = spk_mat;
    np11_struct.cat = [np11_struct.cat spk_mat];
    
end

% 1-shank 2nd session
pts = 1001:5000; %-- UNIQUE TO THIS DATASET; CHANGE BACK LATER (-3000 TO 2000 ms)
np12_struct.cat = []; 
for jj=1:size(np12_pullStarts.unitTimeTrial,3)
    spk_mat = zeros(numel(np12_pullStarts.valid_units),numel(pts));
    
    for mm=1:numel(np12_pullStarts.valid_units)
        if jj==1
            np12_struct.lag(mm) = find(np12_pullStarts.trAVG{np12_pullStarts.valid_units(mm)}(1,pts)==max(np12_pullStarts.trAVG{np12_pullStarts.valid_units(mm)}(1,pts)),1);
            np12_struct.mag(mm) = max(np12_pullStarts.trAVG{np12_pullStarts.valid_units(mm)}(1,pts))-mean(np12_pullStarts.trAVG{np12_pullStarts.valid_units(mm)}(1,pts));
        end
        spk_mat(mm,:) = conv( np12_pullStarts.unitTimeTrial(np12_pullStarts.valid_units(mm),pts,jj) ,spk_kernel,'same');   
    end
    
    np12_struct.wins(jj).win_dat = spk_mat;
    np12_struct.cat = [np12_struct.cat spk_mat];
    
end

% 4-shank np20
pts = 1:4000; %-- UNIQUE TO THIS DATASET; CHANGE BACK LATER (-2000 TO 2000 ms)
np2_struct.cat = []; 
for jj=1:size(np2_pullStarts.unitTimeTrial,3)
    spk_mat = zeros(numel(np2_pullStarts.valid_units),numel(pts));
    
    for mm=1:numel(np2_pullStarts.valid_units)
        if jj==1
            np2_struct.lag(mm) = find(np2_pullStarts.trAVG{np2_pullStarts.valid_units(mm)}(1,pts)==max(np2_pullStarts.trAVG{np2_pullStarts.valid_units(mm)}(1,pts)),1);
            np2_struct.mag(mm) = max(np2_pullStarts.trAVG{np2_pullStarts.valid_units(mm)}(1,pts))-mean(np2_pullStarts.trAVG{np2_pullStarts.valid_units(mm)}(1,pts));
        end
        spk_mat(mm,:) = conv( np2_pullStarts.unitTimeTrial(np2_pullStarts.valid_units(mm),pts,jj) ,spk_kernel,'same');   
    end
    
    np2_struct.wins(jj).win_dat = spk_mat;
    np2_struct.cat = [np2_struct.cat spk_mat];
    
end

% z-score normalization of concatenated data 
for qq=1:size(np11_struct.cat,1)
    np11_struct.cat(qq,:) = (np11_struct.cat(qq,:) - mean(np11_struct.cat(qq,:))) ./ std(np11_struct.cat(qq,:));
end

for qq=1:size(np12_struct.cat,1)
    np12_struct.cat(qq,:) = (np12_struct.cat(qq,:) - mean(np12_struct.cat(qq,:))) ./ std(np12_struct.cat(qq,:));
end

for qq=1:size(np2_struct.cat,1)
    np2_struct.cat(qq,:) = (np2_struct.cat(qq,:) - mean(np2_struct.cat(qq,:))) ./ std(np2_struct.cat(qq,:));
end

% run pca 
[m_np11, mapping_np11] = compute_mapping(np11_struct.cat', 'PCA',3);
[m_np12, mapping_np12] = compute_mapping(np12_struct.cat', 'PCA',3);
[m_np2, mapping_np2] = compute_mapping(np2_struct.cat', 'PCA',3);

%% Plotting
map2 = TNC_CreateRBColormap(2000,'cpb');
map11 = linspace(1,map2(1,1),1000)'; 
map12 = linspace(1,map2(1,2),1000)';
map13 = linspace(1,map2(1,3),1000)';
map1 = [map11, map12, map13]; 
%map1 = TNC_CreateRBColormap(1000,'wblue');

exag_map = [map1; map2]; %[map1(1000:-1:1,:) ; map2];
figure(10); clf;
%subplot(131);
% plot(m_np11(:,1),m_np11(:,2),'k.');
hist3(m_np11(:,1:2),[100 100],'CdataMode','auto','EdgeColor','none'); colormap(exag_map);
%xlabel('PC1')
%ylabel('PC2')
colorbar;
view(2);
caxis([0 3000]);
xlim([-7 15]); ylim([-8 10]);
grid off
set(gca,'tickDir','out')
%set(gca,'xticklabel',{[]})
%set(gca,'yticklabel',{[]})
%title('WR39_100219','Interpreter','none')
print(fullfile('/Users/parkj/Dropbox (HHMI)/np1p0_2p0_compare','pcScore_dSTR_1shank_WR39_100219'),'-dpdf','-bestfit','-painters')

figure(11); clf;
%subplot(132);
% plot(m_np12(:,1),m_np12(:,2),'k.');
hist3(m_np12(:,1:2),[100 100],'CdataMode','auto','EdgeColor','none'); colormap(exag_map);
%xlabel('PC1')
%ylabel('PC2')
%colorbar;
view(2);
caxis([0 3000]);
xlim([-7 15]); ylim([-8 10]);
grid off
set(gca,'tickDir','out')
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
%title('pcScore_dSTR_1shank_WR40_081919','Interpreter','none')
print(fullfile('/Users/parkj/Dropbox (HHMI)/np1p0_2p0_compare','pcScore_dSTR_1shank_WR40_081919'),'-dpdf','-bestfit')

figure(12); clf;
%subplot(133);
% plot(m_np2(:,1),m_np2(:,2),'k.');
hist3(m_np2(:,1:2),[100 100],'CdataMode','auto','EdgeColor','none'); colormap(exag_map);
%xlabel('PC1')
%ylabel('PC2')
%colorbar;
view(2);
caxis([0 3000]); 
xlim([-7 15]); ylim([-8 10]);
grid off
set(gca,'tickDir','out')
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
%title('pcScore_dSTR_4shank_WR39_100819','Interpreter','none')
print(fullfile('/Users/parkj/Dropbox (HHMI)/np1p0_2p0_compare','pcScore_dSTR_4shank_WR39_100819'),'-dpdf','-bestfit')

figure(1); clf;
subplot(131);
scatter(mapping_np11.M(:,1),mapping_np11.M(:,2),30,np11_pullStarts.ypos(np11_pullStarts.valid_units),'filled'); colormap(lat_map);
subplot(132); 
scatter(mapping_np12.M(:,1),mapping_np12.M(:,2),30,np12_pullStarts.ypos(np12_pullStarts.valid_units),'filled'); colormap(lat_map);
subplot(133);
scatter(mapping_np2.M(:,1),mapping_np2.M(:,2),30,np2_pullStarts.ypos(np2_pullStarts.valid_units),'filled'); colormap(lat_map);

[rho_p11,p_np11] = corr([np11_pullStarts.ypos(np11_pullStarts.valid_units)',mapping_np11.M(:,2)]); %m_np11(:,2)]);
[rho_p12,p_np12] = corr([np12_pullStarts.ypos(np12_pullStarts.valid_units)',mapping_np12.M(:,2)]); %m_np11(:,2)]);
[rho_p2,p_np2] = corr([np2_pullStarts.ypos(np2_pullStarts.valid_units)',mapping_np2.M(:,2)]); %m_np2(:,2)]);

[np11_struct.lag_v,np11_struct.lag_i] = sort(np11_struct.lag,'ascend');
[np12_struct.lag_v,np12_struct.lag_i] = sort(np12_struct.lag,'ascend');
[np2_struct.lag_v,np2_struct.lag_i] = sort(np2_struct.lag,'ascend');

[np11_struct.mag_v,np11_struct.mag_i] = sort(np11_struct.mag,'descend');
[np12_struct.mag_v,np12_struct.mag_i] = sort(np12_struct.mag,'descend');
[np2_struct.mag_v,np2_struct.mag_i] = sort(np2_struct.mag,'descend');

[np11_struct.dp_v,np11_struct.dp_i] = sort(np11_pullStarts.ypos(np11_pullStarts.valid_units),'descend');

figure(3); clf; 
imagesc(smooth2a(np11_struct.cat(np11_struct.lag_i,4e3*10:4e3*20),0,10),[-5 5]); colormap(phys_map); hold on;
plot([2e3:4e3:4e3*10],[ones(1,10)],'k*'); box off;
set(gca,'tickDir','out')
print(fullfile('/Users/parkj/Dropbox (HHMI)/np1p0_2p0_compare','zScoreExampleTrace_WR39_100219'),'-dpdf','-bestfit','-painters')

figure(4); clf; 
imagesc(smooth2a(np12_struct.cat(np12_struct.lag_i,4e3*10:4e3*20),0,10),[-5 5]); colormap(phys_map); hold on;
plot([2e3:4e3:4e3*10],[ones(1,10)],'k*'); box off;
set(gca,'tickDir','out')
print(fullfile('/Users/parkj/Dropbox (HHMI)/np1p0_2p0_compare','zScoreExampleTrace_WR40_081919'),'-dpdf','-bestfit','-painters')

figure(5); clf; 
imagesc(smooth2a(np2_struct.cat(np2_struct.lag_i,4e3*10:4e3*20),0,10),[-5 5]); colormap(phys_map); hold on;
plot([2e3:4e3:4e3*10],[ones(1,10)],'k*'); box off;
set(gca,'tickDir','out')
print(fullfile('/Users/parkj/Dropbox (HHMI)/np1p0_2p0_compare','zScoreExampleTrace_WR39_100819'),'-dpdf','-bestfit','-painters')
% exag_map = [lat_map(2000:-1:1,:) ; lat_map(4000:-1:2001,:)];

figure(6); clf;
%subplot(121); 
plt.xpos1 = 375+np11_pullStarts.xpos(np11_pullStarts.valid_units)+randn(1,numel(np11_pullStarts.valid_units))*30;
plt.ypos1 = np11_pullStarts.ypos(np11_pullStarts.valid_units);
plt.seq1 = find(np11_struct.mag>median(np11_struct.mag) & np11_struct.lag<3000);
[np11_struct.lag_v,np11_struct.lag_i] = sort(np11_struct.lag(plt.seq1),'ascend');
plot(plt.xpos1(plt.seq1(np11_struct.lag_i)),plt.ypos1(plt.seq1(np11_struct.lag_i)),'color',[0.85 0.85 0.85]);
hold on;
quiver(plt.xpos1(plt.seq1(np11_struct.lag_i)),plt.ypos1(plt.seq1(np11_struct.lag_i)),[diff(plt.xpos1(plt.seq1(np11_struct.lag_i))) 0],[diff(plt.ypos1(plt.seq1(np11_struct.lag_i))) 0],'color',[0.5 0.5 0.5]);
hold on;
scatter(plt.xpos1,plt.ypos1,2300*np11_struct.mag,np11_struct.lag,'filled'); axis([-100 800 2200 3000]); set(gca,'YDir','reverse'); colormap(exag_map); 
ylim([2200 3100])
xticks([375])
xlim([-100 850])
box off;
print(fullfile('/Users/parkj/Dropbox (HHMI)/np1p0_2p0_compare','activationProbeMap_WR39_100219'),'-dpdf','-bestfit','-painters')

figure(7); clf;
%subplot(121); 
plt.xpos1 = 375+np12_pullStarts.xpos(np12_pullStarts.valid_units)+randn(1,numel(np12_pullStarts.valid_units))*30;
plt.ypos1 = np12_pullStarts.ypos(np12_pullStarts.valid_units);
plt.seq1 = find(np12_struct.mag>median(np12_struct.mag) & np12_struct.lag<3000);
[np12_struct.lag_v,np12_struct.lag_i] = sort(np12_struct.lag(plt.seq1),'ascend');
plot(plt.xpos1(plt.seq1(np12_struct.lag_i)),plt.ypos1(plt.seq1(np12_struct.lag_i)),'color',[0.85 0.85 0.85]);
hold on;
quiver(plt.xpos1(plt.seq1(np12_struct.lag_i)),plt.ypos1(plt.seq1(np12_struct.lag_i)),[diff(plt.xpos1(plt.seq1(np12_struct.lag_i))) 0],[diff(plt.ypos1(plt.seq1(np12_struct.lag_i))) 0],'color',[0.5 0.5 0.5]);
hold on;
scatter(plt.xpos1,plt.ypos1,2300*np12_struct.mag,np12_struct.lag,'filled'); axis([-100 800 2200 3000]); set(gca,'YDir','reverse'); colormap(exag_map); 
ylim([2200 3100])
xticks([375])
xlim([-100 850])
box off;
print(fullfile('/Users/parkj/Dropbox (HHMI)/np1p0_2p0_compare','activationProbeMap_WR40_081919'),'-dpdf','-bestfit','-painters')

figure(8); clf;
%subplot(122); 
plt.xpos = np2_pullStarts.xpos(np2_pullStarts.valid_units)+randn(1,numel(np2_pullStarts.valid_units))*30;
plt.ypos = np2_pullStarts.ypos(np2_pullStarts.valid_units);
plt.seq = find(np2_struct.mag>median(np2_struct.mag) & np2_struct.lag<3000);
[np2_struct.lag_v,np2_struct.lag_i] = sort(np2_struct.lag(plt.seq),'ascend'); 
plot(plt.xpos(plt.seq(np2_struct.lag_i)),plt.ypos(plt.seq(np2_struct.lag_i)),'color',[0.85 0.85 0.85]);
hold on;
quiver(plt.xpos(plt.seq(np2_struct.lag_i)),plt.ypos(plt.seq(np2_struct.lag_i)),[diff(plt.xpos(plt.seq(np2_struct.lag_i))) 0],[diff(plt.ypos(plt.seq(np2_struct.lag_i))) 0],'color',[0.5 0.5 0.5]);
hold on;
scatter(plt.xpos,plt.ypos,3000*np2_struct.mag,np2_struct.lag,'filled'); axis([-100 800 2200 3000]); set(gca,'YDir','reverse'); colormap(exag_map); 
set(gca,'tickDir','out')
xticks([0 250 500 750])
xlim([-100 850])
ylim([2200 3100])
box off;
print(fullfile('/Users/parkj/Dropbox (HHMI)/np1p0_2p0_compare','activationProbeMap_WR39_100819'),'-dpdf','-bestfit','-painters')









