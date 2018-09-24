%% cosine waves (tuning curves)
x = linspace(-pi,pi,100);
base1 = cos(x-pi/2);
base2 = cos(x);
base3 = cos(x+pi/2);

signalVar = 0.05;
noiseAmp = max(abs(cos(x)))*0.1;

for n = 1:30
    if rem(n,3)==1
        for t = 1:150
            dat(t,:,n) = base1+randn(1,length(x))*noiseAmp+randn(1)*signalVar;
        end
    elseif rem(n,3)==2
        for t = 1:150
            dat(t,:,n) = base2+randn(1,length(x))*noiseAmp+randn(1)*signalVar;
        end
    elseif rem(n,3)==0
        for t = 1:150
            dat(t,:,n) = base3+randn(1,length(x))*noiseAmp+randn(1)*signalVar;
        end
    end
end
            
%% run pca
rsDat = reshape(dat,size(dat,1)*size(dat,2),size(dat,3)); 
[loadings, scores] = pca(rsDat); 
pcScoreMat = reshape(scores,150,100,30); 
plot(pcScoreMat(1,:,1)); hold on; plot(pcScoreMat(1,:,2)); plot(pcScoreMat(1,:,3));

%% cosine waves (tuning curves)
x = linspace(-pi,pi,100);
base1 = max(cos(x-pi/2),0);
base2 = max(cos(x+pi/2),0);

signalVar = 0.05;
noiseAmp = max(abs(cos(x)))*0.1;

for n = 1:30
    if rem(n,2)==1
        for t = 1:150
            dat(t,:,n) = base1+randn(1,length(x))*noiseAmp+randn(1)*signalVar;
        end
    elseif rem(n,2)==0
        for t = 1:150
            dat(t,:,n) = base2+randn(1,length(x))*noiseAmp+randn(1)*signalVar;
        end
    end
end
            
%% run pca
rsDat = reshape(dat,size(dat,1)*size(dat,2),size(dat,3)); 
[loadings, scores] = pca(rsDat); 
pcScoreMat = reshape(scores,150,100,30); 
plot(pcScoreMat(1,:,1)); hold on; plot(pcScoreMat(1,:,2)); plot(pcScoreMat(1,:,3));