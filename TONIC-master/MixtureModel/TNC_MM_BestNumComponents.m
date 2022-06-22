%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%
% Purpose: identify the number of components in the Gaussian model 
%          which results in the lowest value of cross-validation criterion
%          among the models with different numbers of Gaussian components

function [ best_mc_cvl ] = TNC_MM_BestNumComponents(ft_file, shank, segment, ...
                                                method, mc_list, ...
                                                num_folds, output_mm_file)
    disp(['ft_file=' ft_file ]);
    disp(['output_mm_file=' output_mm_file ]);
    shank
    method 
    mc_list
    num_folds
    name_prefix = ft_file(1:(length(char(ft_file))-7));
    disp(['mc_list=' num2str(mc_list)]); 
    mc          = str2num(char(strsplit(mc_list, ',')'))';                 
    iShank      = int32(str2num(shank));
    iSeg        = int32(str2num(segment));
    nf          = int32(str2num(char(num_folds)));
    disp(['mc= ' num2str(mc) ' iShank=' num2str(iShank) ' iSeg=' num2str(iSeg) ]);

    mmStruct.program    = 'TNC_MM_BestNumComponents.m';
    mmStruct.method     = method;
    mmStruct.numFolds   = num_folds;
    cvl                 = zeros(1, length(mc)); 
    aic                 = zeros(1, length(mc));
    mml                 = zeros(1, length(mc));
    mml                 = zeros(1, length(mc));
    for m=1:length(mc)      
        input_mm_file = strcat(name_prefix,'_shank', num2str(iShank), ...
                               '_seg', num2str(iSeg), '_mc', num2str(mc(m)), ...
                               '_mm.mat');
        disp(['Processing file ' input_mm_file]);
        mmd = load(input_mm_file);
        cvl1 = mmd.mmStruct.shank(iShank).seg(iSeg).cvl;
        % Cross validation likelihood
        cvl(m) = cvl1;
        disp(['iSeg=' num2str(iSeg) ' iShank=' num2str(iShank) ' mc=' num2str(mc(m)) ' cvl=' num2str(cvl(m))]);
        LL0  = mmd.mmStruct.shank(iShank).seg(iSeg).LL0;
        N    = mmd.mmStruct.shank(iShank).seg(iSeg).N;
        K    = mmd.mmStruct.K;
        Ns   = 2 + K + K*(K+1)/2;
        pi   = 3.142;
        PI   = mmd.mmStruct.shank(iShank).seg(iSeg).PI;
        aic(m) = log(LL0) - mc(m);
        bic(m) = log(LL0) - mc(m)*[log(N) - log(2*pi)]/2.;
        mml(m) = log(LL0) -[mc(m)*(Ns+1) + (mc(m)+Ns)*log(N/12.) + Ns*sum(log(PI))]/2.;
        clear mmd;
    end
    best_m_cvl = min(find(cvl == max(cvl)));
    best_mc_cvl = mc(best_m_cvl);
    disp(['iSeg=' num2str(iSeg) ' best_m_cvl=' num2str(best_m_cvl) ' best_mc_cvl=' num2str(best_mc_cvl)]);

    best_m_aic = min(find(aic == max(aic)));
    best_mc_aic = mc(best_m_aic);
    disp(['iSeg=' num2str(iSeg) ' best_m_aic=' num2str(best_m_aic) ' best_mc_aic=' num2str(best_mc_aic)]);

    best_m_bic = min(find(bic == max(bic)));
    best_mc_bic = mc(best_m_bic);
    disp(['iSeg=' num2str(iSeg) ' best_m_bic=' num2str(best_m_bic) ' best_mc_bic=' num2str(best_mc_bic)]);

    best_m_mml = min(find(mml == max(mml)));
    best_mc_mml = mc(best_m_mml);
    disp(['iSeg=' num2str(iSeg) ' best_m_mml=' num2str(best_m_mml) ' best_mc_mml=' num2str(best_mc_mml)]);

    % Output the best mc to mm file
    mmStruct.shank(iShank).seg(iSeg).best_mc_cvl  = best_mc_cvl;
    mmStruct.shank(iShank).seg(iSeg).best_mc_aic  = best_mc_aic;
    mmStruct.shank(iShank).seg(iSeg).best_mc_bic  = best_mc_bic;
    mmStruct.shank(iShank).seg(iSeg).best_mc_mml  = best_mc_mml;

    mmStruct.shank(iShank).seg(iSeg).mc_range = mc;
    mmStruct.shank(iShank).seg(iSeg).cvl      = cvl;     

    % Load data from the first good representative file
    success = 0;
    fold = 0;
    while success == 0 
        input_mm_file = strcat(name_prefix,'_shank',num2str(iShank), ...
                               '_seg', num2str(iSeg), ...
                               '_mc',  num2str(best_mc_cvl), ...
                               '_mm_', num2str(fold), '.mat');
        try
            mmData = load(input_mm_file);
            assert(numel(mmData.mmStruct.shank(iShank).seg(iSeg).Y) > 0);
            mmStruct.thr                           = mmData.mmStruct.thr;
            mmStruct.K                             = mmData.mmStruct.K;
            mmStruct.shank(iShank).seg(iSeg).N     = N;
            disp(['fold=' num2str(fold) ' N=' num2str(N)]);
            mmStruct.shank(iShank).seg(iSeg).Y_true= mmData.mmStruct.shank(iShank).seg(iSeg).Y_true;
            mmStruct.shank(iShank).seg(iSeg).Y     = mmData.mmStruct.shank(iShank).seg(iSeg).Y;
            mmStruct.shank(iShank).seg(iSeg).Q     = mmData.mmStruct.shank(iShank).seg(iSeg).Q;
            mmStruct.shank(iShank).seg(iSeg).classError = mmData.mmStruct.shank(iShank).seg(iSeg).classError;
            mmStruct.shank(iShank).seg(iSeg).LL    = mmData.mmStruct.shank(iShank).seg(iSeg).LL;
            mmStruct.shank(iShank).seg(iSeg).UC    = mmData.mmStruct.shank(iShank).seg(iSeg).UC;
            mmStruct.shank(iShank).seg(iSeg).PI    = mmData.mmStruct.shank(iShank).seg(iSeg).PI;
            mmStruct.shank(iShank).seg(iSeg).MU    = mmData.mmStruct.shank(iShank).seg(iSeg).MU;
            mmStruct.shank(iShank).seg(iSeg).SIGMA = mmData.mmStruct.shank(iShank).seg(iSeg).SIGMA;
            mmStruct.shank(iShank).seg(iSeg).TSM   = mmData.mmStruct.shank(iShank).seg(iSeg).TSM;
            success = 1;
        catch
            disp(['file ' input_mm_file ' is not good']);
        end
        fold = fold + 1;
        if fold > double(num_folds)
            break;
        end
    end

    % Save results
    disp(['save ' output_mm_file ' mmStruct']);
    eval(['save ' output_mm_file ' mmStruct']); 
