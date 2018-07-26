%% TRAIN A decoder of movement from neural data

train=0;
%(should pick leave out trials for cross validation now and then check against that subset later)


train_i = spk.dens_dec(1).dens_rew;
test_i  = spk.dens_dec(2).dens_rew;
train_o = [sgolayfilt( spk.dens_dec(1).mv_raw-mean(spk.dens_dec(1).mv_raw(1:500)) , 3 , 101 ) ; spk.lk_raw.*1e3];
test_o  = [sgolayfilt( spk.dens_dec(2).mv_raw-mean(spk.dens_dec(2).mv_raw(1:500)) , 3 , 101 ) ; spk.lk_raw.*1e3];
xsamps = -2000:8000;
figure(21); clf;
subplot(221); shadedErrorBar(xsamps,train_o(1,:),spk.mv_raw_E,{'color',[0  0.33 0.5]}); hold on; title('Training data');
subplot(223); shadedErrorBar(xsamps,test_o(1,:),spk.mv_raw_E,{'color',[0  0.33 0.5]}); hold on; title('Hold out');
subplot(222); shadedErrorBar(xsamps,train_o(2,:),spk.lk_raw_E.*1e3,{'color',[0.5 0 0]}); hold on; title('Training data');
subplot(224); shadedErrorBar(xsamps,train_o(2,:),spk.lk_raw_E.*1e3,{'color',[0.5 0 0]}); hold on; title('Hold out');
drawnow;

if train
   nn      = [33 6 4 2];
   dIn     = [1];
   dIntern = [1];
   dOut    = [1];
   net     = CreateNN(nn,dIn,dIntern,dOut);
   netLM2 = train_LM( train_i , train_o , net, 25, 15);

   decode_train = NNOut(train_i,netLM2);
   decode_test = NNOut(test_i,netLM2);
end

figure(21);
subplot(221); plot(xsamps,decode_train(1,:),'linewidth',2,'color',[0 0.66 1]);
subplot(223); plot(xsamps,decode_test(1,:),'linewidth',2,'color',[0 0.66 1]);
subplot(222); plot(xsamps,decode_train(2,:),'linewidth',2,'color',[1 0 0]);
subplot(224); plot(xsamps,decode_test(2,:),'linewidth',2,'color',[1 0 0]);