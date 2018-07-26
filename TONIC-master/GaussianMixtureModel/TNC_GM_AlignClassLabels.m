%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

% This function performs global alignment of 2 strings                  

function [classError] = TNC_GM_AlignClassLabels(String1,String2,verbose)
    Match=1;
    Sub  =-1;
    Gap  =-2;

    Size1 = length(String1)+1;
    Size2 = length(String2)+1;

    SimMatrix = populate_similarity_matrix(Size1,Size2,String1,String2,...
                                           Match,Sub,Gap);
    Alignment = compute_best_alignment(Size1,Size2,String1,String2,SimMatrix,...
                                            Match,Sub,Gap);
    classError.tot = sum(Alignment.AlignSymbols  == ' ')/length(String2);
    classError.ins   = sum(Alignment.GappedString2 == '-')/length(String2); 
    classError.del   = sum(Alignment.GappedString1 == '-')/length(String2);
    classError.sub   = classError.tot - classError.ins - classError.del;

    if verbose
        output_alighment(Alignment)
    end

%-------------------------------------------------------------------------------
    
function SimMatrix = populate_similarity_matrix(Size1,Size2,String1,String2,...
                                                Match,Sub,Gap)
    SimMatrix = zeros(Size1, Size2); 
    for i=2:Size1
        SimMatrix(i,1) = SimMatrix(i-1,1)+Gap;
    end
    for j=2:Size2
        SimMatrix(1,j) = SimMatrix(1,j-1)+Gap;
    end 
    k = 2;
    while k <= Size1 && k <= Size2
        if k <= Size1
            for i=k:Size1
                if String1(i-1)==String2(k-1)
                    p = Match;
                else
                    p = Sub;
                end
                SimMatrix(i,k) = max(max(SimMatrix(i  ,k-1)+Gap,...
                                         SimMatrix(i-1,k-1)+p), ...
                                         SimMatrix(i-1,k  )+Gap);
            end
        end
        if k <= Size2 
            for j=k:Size2
                if String1(k-1)==String2(j-1)
                    p = Match;
                else
                    p = Sub;
                end
                SimMatrix(k,j) = max(max(SimMatrix(k-1,j  )+Gap,...
                                         SimMatrix(k-1,j-1)+p), ...
                                         SimMatrix(k  ,j-1)+Gap);
            end
        end
        k = k+1;
    end

%-------------------------------------------------------------------------------

function Alignment = compute_best_alignment(Size1,Size2,String1,String2,SimMatrix,...
                                            Match,Sub,Gap)
    Alignment.GappedString1 = [];
    Alignment.GappedString2 = [];
    Alignment.AlignSymbols  = [];
    i = Size1;
    j = Size2;
    while i > 1 || j > 1
        if i > 1 && j > 1 && String1(i-1) == String2(j-1)
            p = Match;
        else
            p = Sub;
        end
        if i > 1 && SimMatrix(i,j) == SimMatrix(i-1,j) + Gap
            Alignment.GappedString1 = [String1(i-1) Alignment.GappedString1];
            Alignment.GappedString2 = ['-'          Alignment.GappedString2];
            Alignment.AlignSymbols  = [' '          Alignment.AlignSymbols];
            i = i-1;
        elseif i > 1 && j > 1 && SimMatrix(i,j) == SimMatrix(i-1,j-1) + p
            Alignment.GappedString1 = [String1(i-1) Alignment.GappedString1];
            Alignment.GappedString2 = [String2(j-1) Alignment.GappedString2];
            if p == Match
                Alignment.AlignSymbols  = ['|' Alignment.AlignSymbols];
            else
                Alignment.AlignSymbols  = [' ' Alignment.AlignSymbols];
            end
            i = i-1;
            j = j-1;
        else
            Alignment.GappedString1 = ['-'          Alignment.GappedString1];
            Alignment.GappedString2 = [String2(j-1) Alignment.GappedString2];
            Alignment.AlignSymbols  = [' '          Alignment.AlignSymbols];
            j = j-1;
        end
    end

%-------------------------------------------------------------------------------

function [] = output_alighment(Alignment)
    disp(Alignment.GappedString1);
    disp(Alignment.AlignSymbols);
    disp(Alignment.GappedString2);
 
