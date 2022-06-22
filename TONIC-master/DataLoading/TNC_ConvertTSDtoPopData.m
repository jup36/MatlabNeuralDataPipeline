function [PopData] = TNC_ConvertTSDtoPopData(source,sessNum)

ucnt = 0;
figure(1); clf;

d = dir([source '*_tsd.mat']);

if numel(d)>0

    for i=1:numel(d)

        disp(['Loading data from >>> ' d(i).name]);
        s = load(d(i).name);

        numUnits = numel(s.shank.unit);

        for k=1:numUnits
            
            ucnt = ucnt+1;
            
            nonZeros = find(diff(s.shank.unit(k).ts)>1);            
            
            PopData.session(sessNum).unit(ucnt).inds    = s.shank.unit(k).inds(nonZeros);            
            PopData.session(sessNum).unit(ucnt).ts      = round(PopData.session(sessNum).unit(ucnt).inds./30);
%             PopData.session(sessNum).unit(ucnt).ts      = s.shank.unit(k).ts(nonZeros);
            PopData.session(sessNum).unit(ucnt).sh      = str2num(d(i).name(strfind(d(i).name,'shank')+5));
            
            [isi] = TNC_QuantISI(PopData.session(sessNum).unit(ucnt).ts);
            PopData.session(sessNum).unit(ucnt).isi     = isi;

        end

    end

    disp(['Session number ' num2str(sessNum) ' was written to the PopData structure and contains ' num2str(ucnt) ' single units.']);    

    dims = ceil(sqrt(ucnt));
        
    for m=1:ucnt            
        figure(1);
%         subplot(dims,dims,m);
        subplot(1,ucnt,m); hold off;
        PopData.isiDist(m,:) = PopData.session(sessNum).unit(m).isi.hist.logCount;
        semilogx(PopData.session(sessNum).unit(m).isi.hist.logTimes,PopData.session(sessNum).unit(m).isi.hist.logCount,'LineWidth',1,'Color',[(m./ucnt) 0.5 1-(m./ucnt)]); hold on;
        semilogx([1 1],[0 max(PopData.session(sessNum).unit(m).isi.hist.logCount)],'k--');
        title(['sh' num2str(PopData.session(sessNum).unit(m).sh) 'u' num2str(m) '   |   ' num2str(round(1000./PopData.session(sessNum).unit(m).isi.stats.mean)) ' Hz ']);
        grid on; axis([1e-2 1e4 0 max(PopData.session(sessNum).unit(m).isi.hist.logCount)]); axis off;
    end        
    
    [ids] = kmeans(PopData.isiDist,3);
    [ys,is] = sort(ids);
    PopData.isiClass.y = ys;
    PopData.isiClass.i = is;
    
    [cMap] = TNC_CreateRBColormap(12,'mbr');
    figure(2); 
    subplot(211); imagesc(PopData.isiDist(is,:)); colormap(cMap);
    subplot(212); imagesc(corr(PopData.isiDist(is,:)'),[0 1]); colormap(cMap);
    
else

    disp('Could not find any sorted files with that name.');

end

