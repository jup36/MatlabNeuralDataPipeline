%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%
% Purpose: identify the number of components in the Gaussian model 
%          which results in the lowest value of cross-validation criterion
%          among the models with different numbers of Gaussian components

function [ cvl ] = TNC_MM_CrossValidation(ft_file, shank, segment_list, ...
                                          method, mc, num_folds, ...
                                          output_mm_file)
    name_prefix = ft_file(1:(length(char(ft_file))-7));
    seg_ids     = str2num(char(strsplit(segment_list,  ',')'))';
    iShank      = int32(str2num(shank));
    nf          = int32(str2num(num_folds));
    disp(['mc= ' num2str(mc)]);
    disp(['seg_ids=' num2str(seg_ids)]);

    mmStruct.program  = 'TNC_MM_CrossValidation.m';
    mmStruct.method   = method;
    mmStruct.mc       = mc;
    mmStruct.numFolds = nf;
    cvl                 = zeros(1, length(seg_ids)); 
    for i=1:length(seg_ids)
        iSeg = int32(seg_ids(i));
        num_good_folds = 0;
        for f =0:nf
            input_mm_file = strcat(name_prefix,'_shank', num2str(iShank), ...
                                               '_seg',   num2str(iSeg), ...
                                               '_mc',    num2str(mc), ...
                                               '_mm_',   num2str(f), '.mat');
            file_exists = (exist(input_mm_file, 'file') == 2);
            if file_exists
                disp(['Processing file ' input_mm_file ' file_exists=' num2str(file_exists)]);
                mmd = load(input_mm_file);
                if f == 0
                    mmStruct.shank(iShank).seg(iSeg).LL0       = mmd.mmStruct.shank(iShank).seg(iSeg).LL;
                else
                    mmStruct.thr = mmd.mmStruct.thr;
                    mmStruct.K   = mmd.mmStruct.K;
                    LLf                                        = mmd.mmStruct.shank(iShank).seg(iSeg).LLf;
                    mmStruct.shank(iShank).seg(iSeg).LLf       = LLf;
                    mmStruct.shank(iShank).seg(iSeg).N         = mmd.mmStruct.shank(iShank).seg(iSeg).N;
                    mmStruct.shank(iShank).seg(iSeg).PI        = mmd.mmStruct.shank(iShank).seg(iSeg).PI;
                    mmStruct.shank(iShank).seg(iSeg).MU        = mmd.mmStruct.shank(iShank).seg(iSeg).MU;
                    mmStruct.shank(iShank).seg(iSeg).SIGMA     = mmd.mmStruct.shank(iShank).seg(iSeg).SIGMA;
                    mmStruct.shank(iShank).seg(iSeg).TSM       = mmd.mmStruct.shank(iShank).seg(iSeg).TSM;
                    mmStruct.shank(iShank).seg(iSeg).UC        = mmd.mmStruct.shank(iShank).seg(iSeg).UC;
                    mmStruct.shank(iShank).seg(iSeg).Y_true    = mmd.mmStruct.shank(iShank).seg(iSeg).Y_true;
                    mmStruct.shank(iShank).seg(iSeg).Y         = mmd.mmStruct.shank(iShank).seg(iSeg).Y;
                    mmStruct.shank(iShank).seg(iSeg).classError= mmd.mmStruct.shank(iShank).seg(iSeg).classError;
                    mmStruct.shank(iShank).seg(iSeg).LL        = mmd.mmStruct.shank(iShank).seg(iSeg).LL;
                    mmStruct.method                            = mmd.mmStruct.method;
                    mmStruct.thr                               = mmd.mmStruct.thr; 
                    disp(['iSeg=' num2str(iSeg) ' iShank=' num2str(iShank) ...
                          ' mc=' num2str(mc) ' iFold=' num2str(f) ...
                          ' LLf=' num2str(mmd.mmStruct.shank(iShank).seg(iSeg).LLf)]);
                    if imag(LLf) == 0 && LLf > -Inf
                        cvl(i) = cvl(i) + LLf;
                        if LLf ~= 0
                            num_good_folds = num_good_folds + 1;
                        end
                    else 
                        break;
                        cvl(i) = -Inf;
                    end
                end
                clear mmd;
            else
                mmStruct.shank(iShank).seg(iSeg).LL0       = -Inf;
                disp(['Warning: input file ' input_mm_file ' does not exist ']);
            end
        end
        if num_good_folds > 0 
            cvl(i) = cvl(i)/num_good_folds;
        end
        if cvl(i) == 0
            cvl(i) = -Inf;
        end
        mmStruct.shank(iShank).seg(iSeg).cvl = cvl(i);
    end

    % Output the best mc to mm file
    mmStruct.shank(iShank).cvl = cvl;        

    % Save results
    eval(['save ' output_mm_file ' mmStruct']); 
    disp(['save ' output_mm_file ' mmStruct']);
