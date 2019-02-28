
%% Case COBEM 2019
% Comparative study of reservoir zoning through FZI and FZI*
% for UNISIM 1, SPE and NORNE models

% class instantiation
d = DirManager(); 

% data directory
dataDir = d.createTempDataDir(mfilename);   

%% Load G and PROPS from chosen model

model = 'spe';

switch model 
    
    case 'unisim1'                
        [G,~] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');
        
                    
    case 'spe'        
        load('../benchmarks/spe10/mrst/G-SPE10.mat');
        load('../benchmarks/spe10/mrst/rock-SPE10.mat');
        PHI = rock.poro; 
        KX  = rock.perm(:,1);
        KY  = rock.perm(:,2);
        KZ  = rock.perm(:,3);        
                
    case 'norne'
        % For NORNE case, we need to recover the original 'grdecl' 
        % (ECLIPSE processed deck file) and fill PROPS with the 
        % original vectors to be inversely mapped later.
        
        load('../benchmarks/norne/mrst/G-NORNE.mat');
        load('../benchmarks/norne/mrst/rock-NORNE.mat');
        load('../benchmarks/norne/mrst/grdecl-NORNE.mat');
        PROPS.PHI = grdecl.PORO; 
        PROPS.KX  = grdecl.PERMX;
        PROPS.KY  = grdecl.PERMY;
        PROPS.KZ  = grdecl.PERMZ;

end

%% Compute parameters 

PHIZ = PHI./(1.0 - PHI);              % normalized porosity
KN = sqrt(KX.^2 + KY.^2 + KZ.^2);     % normalized permeability

% reservoir quality index (RQI)
% needs 'isinf' to remove null porosity cells
% The conversion factor 0.0314 is not necessary for 
% processed data by MRST because they are already in milidarcy

r = @(var) sqrt(var./PHI);
RQIN = r(KN);  RQIN(isinf(RQIN)) = 0; % normalized

% flow zone indicator (FZI)
FZIN = RQIN./PHIZ; FZIN(isnan(FZIN)) = 0; % normalized

% flow zone indicator * (FZI*)
%?from Paiaman (2015), DOI:10.1002/ente.201500010
FZINStar = RQIN;


% discrete rock type (DRT) - neperian log
r = @(FZI) round(2*log(FZI) + 10.6);
DRTN_LN = r(FZIN); DRTN_LN(isinf(DRTN_LN)) = 0; % normalized

% DRT*, for FZI*
DRTNStar_LN = r(FZINStar); DRTNStar_LN(isinf(DRTNStar_LN)) = 0; % normalized



% P = computeParams(G,PROPS);
% 
% field statistics for chosen properties
% SN = printStats(P,{'DRTN_LN','DRTNStar_LN'},'n');

