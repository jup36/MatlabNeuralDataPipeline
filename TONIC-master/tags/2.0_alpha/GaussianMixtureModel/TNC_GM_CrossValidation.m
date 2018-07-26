%
% Copyright (C) 2013 by Howard Hughes Medical Institute.
%
% Purpose: identify the number of components in the Gaussian model 
%          which results in the lowest value of cross-validation criterion
%          among the models with different numbers of Gaussian components

function [ cvl ] = TNC_GM_CrossValidation(ft_file, shank, segment_list, ...
                                          method, model, mc, num_folds, ...
                                          output_gm_file)
    name_prefix = ft_file(1:(length(char(ft_file))-7));
    seg_ids     = str2num(char(strsplit(segment_list,  ',')'))';
    iShank      = int32(str2num(shank));
    nf          = int32(str2num(num_folds));
    disp(['mc= ' num2str(mc)]);
    disp(['seg_ids=' num2str(seg_ids)]);

    gmStruct.program    = 'TNC_GM_CrossValidation.m';
    gmStruct.method     = method;
    gmStruct.model      = model;
    gmStruct.mc         = mc;
    cvl                 = zeros(1, length(seg_ids)); 
    for i=1:length(seg_ids)
        iSeg = int32(seg_ids(i));
        num_good_folds = 0;
        for f =1:nf
            input_gm_file = strcat(name_prefix,'_shank', num2str(iShank), ...
                                               '_seg',   num2str(iSeg), ...
                                               '_mc',    num2str(mc), ...
                                               '_gm_',   num2str(f), '.mat');
            disp(['Processing file ' input_gm_file]);
            gmd = load(input_gm_file);
            gmStruct.snr = gmd.gmStruct.snr;
            LLf = gmd.gmStruct.shank(iShank).seg(iSeg).LLf;
            gmStruct.shank(iShank).seg(iSeg).LLf = LLf;
            disp(['iSeg=' num2str(iSeg) ' iShank=' num2str(iShank) ' mc=' num2str(mc) ' iFold=' num2str(f) ' LLf=' num2str(LLf)]);
            if imag(LLf) == 0 && LLf > -Inf
                cvl(i) = cvl(i) + LLf;
                if LLf ~= 0
                    num_good_folds = num_good_folds + 1;
                end
            else 
                break;
                cvl(i) = -Inf;
            end
            clear gmd;
        end
        if num_good_folds > 0 && cvl(i) > -Inf && cvl(i) ~= 0
            cvl(i) = cvl(i)/num_good_folds;
        end
        gmStruct.shank(iShank).seg(iSeg).cvl = cvl(i);
    end

    % Output the best mc to gm file
    gmStruct.shank(iShank).cvl = cvl;        

    % Save results
    eval(['save ' output_gm_file ' gmStruct']); 
    disp(['save ' output_gm_file ' gmStruct']);
