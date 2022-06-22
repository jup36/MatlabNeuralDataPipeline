function [lfp] = TNC_ExtractLFP(fileName,chan_list)

    %=============================
    % USE ELLIPTICAL FILTER
    par.sr = 30000;
    par.detect_fmin = 0.5;
    par.detect_fmax = 200;


    %=============================
    %=============================
    disp(['0) Loading the raw data...']);
%     eList = 'e: ';
%     for i=1:numel(chan_list)
%         eList = [eList num2str(chan_list(i))];
%         if i<numel(chan_list)
%             eList = [eList ',']        
%         end
%     end
    Ns5DATA = openNSx('report','read',fileName);
    
    %=============================
    %=============================
    disp(['1) Filtering the raw data...']);
    lfpData = zeros( numel(chan_list) , size(Ns5DATA.Data,2) , 'double');
    for i=1:numel(chan_list)
        [lowBandData,hiBandData]    = TNC_FilterData(double(Ns5DATA.Data(chan_list(i),:)),30000,0,0,[0 1]);
        lfpData(i,:)                = lowBandData.values;
    end    
    lfp.lowCutOff   = lowBandData.lowCutOff;
    lfp.hiCutOff    = lowBandData.hiCutOff;
    lfp.rawSampleRate  = lowBandData.sampleRate;

    %=============================
    %=============================
    disp(['2) Common mode rejection, downsampling...']);
    medVals         = median(lfpData,1);
    lfp.trode(8).sh = 0;
    
    decimate_ratio  = 30;
    for i=1:numel(chan_list)
        lfp.trode(i).v     = decimate(lfpData(i,:)-medVals , decimate_ratio) ;
        lfp.trode(i).vz    = ( lfp.trode(i).v - mean(lfp.trode(i).v) ) ./ std(lfp.trode(i).v);
        lfp.trode(i).ch    = chan_list(i);
        lfp.trode(i).sh    = i;
    end
    
    lfp.sampleRate = lfp.rawSampleRate ./ decimate_ratio;

