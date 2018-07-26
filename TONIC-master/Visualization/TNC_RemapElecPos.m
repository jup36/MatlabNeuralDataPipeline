function [row,col,rSp,cSp,electrodeMatrix] = TNC_RemapElecPos(electrodeNumber,arrayType)
% FUNCTION DETAILS: For a given style of array this function will remap electrode numbers into row and column position on the silicon probe arrray. Together with row spacing and column spacing this can be used to create physical maps of the electrode arrays or associated activity on the arrays
% _________________________________________________________________________
% PART OF THE TONIC PACKAGE
%   developed by JOSHUA DUDMAN
%   begun at COLUMBIA UNIVERSITY, continuing at HHMI / JFRC
% 
% BUG REPORTING: dudman.web@gmail.com
% CONTRIBUTIONS: dudmanlab.org
% _________________________________________________________________________
%
% TO DO: standard naming of array types

switch arrayType
    
    case 'NN_w64'
        electrodeMatrix = [ 28 20 12  5 61 58 50 42 ;
                            23 15  7  4 60 53 45 37 ;
                            21 13  3 30 63 51 43 36 ;
                            29 22 14  2 35 62 52 44 ;
                            26 18 10 31 59 56 48 40 ;
                            25 17  9  6 34 55 47 39 ;
                            19 11  1 32 57 49 41 38 ;
                            27 24 16  8 33 64 54 46 ];
        
        [row,col] = find(electrodeMatrix==electrodeNumber);
        rSp = 200;
        cSp = 200;
        
    case 'NN_b64'
        electrodeMatrix = [ 27 24 16  8 33 64 54 46 ;
                            19 11  1 32 57 49 41 38 ;
                            25 17  9  6 34 55 47 39 ;
                            26 18 10 31 59 56 48 40 ;
                            29 22 14  2 35 62 52 44 ;
                            21 13  3 30 63 51 43 36 ;
                            23 15  7  4 60 53 45 37 ;
                            28 20 12  5 61 58 50 42 ];
        [row,col] = find(electrodeMatrix==electrodeNumber);
        rSp = 20;
        cSp = 200;

        
    case 'NN_b32'
        electrodeMatrix = [ 27 24 16  8  ;
                            19 11  1 32  ;
                            25 17  9  6  ;
                            26 18 10 31  ;
                            29 22 14  2  ;
                            21 13  3 30  ;
                            23 15  7  4  ;
                            28 20 12  5  ];
        [row,col] = find(electrodeMatrix==electrodeNumber);
        rSp = 20;
        cSp = 200;        
        
    case 'NN_b64_dhs'        
        electrodeMatrix = [ [1:8]' 8+[1:8]' 16+[1:8]' 24+[1:8]' 32+[1:8]' 40+[1:8]' 48+[1:8]' 56+[1:8]' ];
        [row,col] = find(electrodeMatrix==electrodeNumber);
        rSp = 20;
        cSp = 200;
        
    case 'NN_b32_dhs'        
        electrodeMatrix = [ 21	28	16	4;
                            29	17	5	12;
                            23	31	14	7;
                            26	19	2	10;
                            22	32	15	6;
                            27	20	3	11;
                            25	30	13	9;
                            24	18	1	8 ];
        [row,col] = find(electrodeMatrix==electrodeNumber);
        rSp = 20;
        cSp = 200;
        
	case 'NN_h64_dhs'        
			electrodeMatrix = [ [1:16]'  16+[1:16]' 32+[1:16]' 48+[1:16]'];
			[row,col] = find(electrodeMatrix==electrodeNumber);
			rSp = 20;
			cSp = 200;
                
    case 'NN_w32'
        disp('Not yet implemented');

    case 'KR_w32'
        disp('Not yet implemented');

    case 'AP_d81'
        [row,col] = ind2sub([9 9], electrodeNumber);
        rSp = 10;
        cSp = 20;
        electrodeMatrix = [];
        
    case 'tetrode'
        electrodeMatrix = [1 ; 2 ; 3 ; 4];
        [row,col] = find(electrodeMatrix==electrodeNumber);
        rSp = 70;   % default value for NeuroCube
        cSp = 0;    % does not matter
        
    otherwise
        row = 1;
        col = electrodeNumber;
        rSp = 200;
        cSp = 200;
        electrodeMatrix = [1:256];

end
