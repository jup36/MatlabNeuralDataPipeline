% plotting sandbox

[mapName] = TNC_CreateRBColormap(1024);

numTrials = size(tonicDataStructure.events.US.ts,1); 

previousTrialInterval(1)    = 0;
previousTrialInterval(2:numTrials) = tonicDataStructure.events.US.ts(1:numTrials-1)-tonicDataStructure.events.CS.ts(1:numTrials-1);

[spkCntPerTrialCS] = TNC_QuantSpksPerTrial(tonicDataStructure.unit(3).responseCSS.image.aligned,[1e4,1.02e4]);
[spkCntPerTrialUS] = TNC_QuantSpksPerTrial(tonicDataStructure.unit(3).responseUSS.image.aligned,[1e4,1.02e4]);

covMat      = cov([spkCntPerTrialCS,spkCntPerTrialUS]');
idx         = kmeans(covMat,2);
sortInds    = [find(idx==1);find(idx==2);find(idx==3)];
covMatSort  = covMat(sortInds,sortInds);

figure(2);
colormap(mapName);
subplot(151);   imagesc(tonicDataStructure.events.EL.responseCSS.image.aligned(sortInds,:));
title('Every Lick aligned to CS');
axis([0.9e4 1.5e4 1 numTrials]);
subplot(152);   imagesc(tonicDataStructure.unit(3).responseCSS.image.aligned(sortInds,:));
title('Unit3 aligned to CS');
axis([0.9e4 1.5e4 1 numTrials]);
subplot(153);   imagesc(tonicDataStructure.unit(3).responseUSS.image.aligned(sortInds,:));
title('Unit3 aligned to US');
axis([0.9e4 1.5e4 1 numTrials]);
subplot(154);   imagesc(tonicDataStructure.unit(1).responseCSS.image.aligned(sortInds,:));
title('Unit1 aligned to CS');
axis([0.9e4 1.5e4 1 numTrials]);
subplot(155);   imagesc(tonicDataStructure.unit(1).responseUSS.image.aligned(sortInds,:));
title('Unit1 aligned to US');
axis([0.9e4 1.5e4 1 numTrials]);

figure(6);
colormap(mapName);
subplot(4,1,1:3);
imagesc(covMatSort);
subplot(414);
plot(1:numTrials,idx,'ro',1:numTrials,(previousTrialInterval-500)./1000,'k.-');
axis([1 numTrials 0 4]);
