%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%
% Purpose: identify the number of components in the Gaussian model 
%          which results in the lowest value of cross-validation criterion
%          among the models with different numbers of Gaussian components

function [ best_mc ] = TNC_GM_BestNumComponents(ft_file, shank, segment, ...
                                                method, model, mc_list, ...
                                                num_folds, output_gm_file)
    disp(['ft_file=' ft_file ]);
    disp(['output_gm_file=' output_gm_file ]);
    shank
    method 
    model
    mc_list
    num_folds
    name_prefix = ft_file(1:(length(char(ft_file))-7));
    disp(['mc_list=' num2str(mc_list)]); 
    mc          = str2num(char(strsplit(mc_list, ',')'))';                 
    iShank      = int32(str2num(shank));
    iSeg        = int32(str2num(segment));
    nf          = int32(str2num(char(num_folds)));
    disp(['mc= ' num2str(mc) ' iShank=' num2str(iShank) ' iSeg=' num2str(iSeg) ]);

    gmStruct.program    = 'TNC_GM_BestNumComponents.m';
    gmStruct.method     = method;
    gmStruct.model      = model;
    gmStruct.numFolds   = num_folds;
    cvl                 = zeros(1, length(mc)); 
    for m=1:length(mc)      
        input_gm_file = strcat(name_prefix,'_shank', num2str(iShank), ...
                               '_seg', num2str(iSeg), '_mc', num2str(mc(m)), ...
                               '_gm.mat');
        disp(['Processing file ' input_gm_file]);
        gmd = load(input_gm_file);
        cvl1 = gmd.gmStruct.shank(iShank).seg(iSeg).cvl;
        cvl(m) = cvl1;
        disp(['iSeg=' num2str(iSeg) ' iShank=' num2str(iShank) ' mc=' num2str(mc(m)) ' cvl=' num2str(cvl(m))]);
        clear gmd;
    end
    best_m = min(find(cvl == max(cvl)));
    disp(['best_m=' num2str(best_m)]);
    best_mc = mc(best_m);
    disp(['iSeg=' num2str(iSeg) ' best_mc=' num2str(best_mc)]);

    % Output the best mc to gm file
    gmStruct.shank(iShank).seg(iSeg).best_mc  = best_mc;
    gmStruct.shank(iShank).seg(iSeg).mc_range = mc;
    gmStruct.shank(iShank).seg(iSeg).cvl      = cvl;     

    input_gm_file0 = strcat(name_prefix,'_shank',num2str(iShank), ...
                                        '_seg',  num2str(iSeg), ...   
                                        '_mc',   num2str(best_mc),'_gm_0.mat');
    gmData0 = load(input_gm_file0);
    try
        gmStruct.snr                           = gmData0.gmStruct.snr;
    catch
        disp('gmData0.gmStruct.snr not available');
    end
    gmStruct.shank(iShank).seg(iSeg).Y_true= gmData0.gmStruct.shank(iShank).seg(iSeg).Y_true;
    gmStruct.shank(iShank).seg(iSeg).Y     = gmData0.gmStruct.shank(iShank).seg(iSeg).Y;
    gmStruct.shank(iShank).seg(iSeg).Q     = gmData0.gmStruct.shank(iShank).seg(iSeg).Q;
    gmStruct.shank(iShank).seg(iSeg).classError = gmData0.gmStruct.shank(iShank).seg(iSeg).classError;
    gmStruct.shank(iShank).seg(iSeg).LL    = gmData0.gmStruct.shank(iShank).seg(iSeg).LL;
    gmStruct.shank(iShank).seg(iSeg).xpar  = gmData0.gmStruct.shank(iShank).seg(iSeg).xpar;
    gmStruct.shank(iShank).seg(iSeg).PI    = gmData0.gmStruct.shank(iShank).seg(iSeg).PI;
    gmStruct.shank(iShank).seg(iSeg).MU    = gmData0.gmStruct.shank(iShank).seg(iSeg).MU;
    gmStruct.shank(iShank).seg(iSeg).SIGMA = gmData0.gmStruct.shank(iShank).seg(iSeg).SIGMA;
    gmStruct.shank(iShank).seg(iSeg).SIGMA_INV ...
                                           = gmData0.gmStruct.shank(iShank).seg(iSeg).SIGMA_INV;
    if gmData0.gmStruct.model == 2
        gmData0.gmStruct.shank(iShank).seg(iSeg).EPS;
    end

    % Determine the best M across all segments
    best_m  = min(find(cvl == max(cvl)));
    best_mc = mc(best_m);
    disp(['All segments: best_mc=' num2str(best_mc)]);
    gmStruct.shank(iShank).best_mc  = best_mc;
    gmStruct.shank(iShank).mc_range = mc;

    % Save results
    disp(['save ' output_gm_file ' gmStruct']);
    eval(['save ' output_gm_file ' gmStruct']); 
