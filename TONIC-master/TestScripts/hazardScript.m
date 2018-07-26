% use this script to calculate hazard functions for t5, t500 and t1000 in
% our task...

t50.pdf = TNC_CreateGaussian(2500,50,5000,1);
t50.cdf = cumsum(t50.pdf);

t500.pdf = TNC_CreateGaussian(2500,250,5000,1);
t500.cdf = cumsum(t500.pdf);

t1000.pdf = TNC_CreateGaussian(2500,500,5000,1);
t1000.cdf = cumsum(t1000.pdf);

size(t50.pdf)
size(t50.cdf)

hazard(1,:) = t50.pdf./(1-t50.cdf);
hazard(2,:) = t500.pdf./(1-t500.cdf);
hazard(3,:) = t1000.pdf./(1-t1000.cdf);


figure(10);
subplot(411);
plot(1:5000,t50.pdf./max(t50.pdf),'k',1:5000,t50.cdf,'r');
subplot(412);
plot(1:5000,t500.pdf./max(t500.pdf),'k',1:5000,t500.cdf,'r');
subplot(413);
plot(1:5000,t1000.pdf./max(t1000.pdf),'k',1:5000,t1000.cdf,'r');
subplot(414);
plot(1:1:4500,hazard(1,1:1:4500)./max(hazard(1,1:1:4500)),'r',1:1:4500,hazard(2,1:1:4500)./max(hazard(2,1:1:4500)),'b',1:1:4500,hazard(3,1:1:4500)./max(hazard(3,1:1:4500)),'k');