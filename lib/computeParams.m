function P = computeParams(G,PROPS)
%COMPUTEPARAMS Compute petrophysical and analytical parameters required to
%              reservoir modelling

% SYNOPSIS: 
%   P = computeParams(G,PROPS);
%
% PARAMETERS:
%   G       - reservoir grid
%   PROPS   - reservoir properties
%
% RETURNS:
%  P     - structure having petrophysical and analytical parameters
%          as 3D matrices (required to keep Cartesian cell indexing)

%{
    Developed at LaMEP/UFPB, Brazil
    @gcpeixoto

%}

% checking
p = {'PHI','KX','KY','KZ'};
if ~isstruct(PROPS) 
    error('wellps:computeParams','PROPS is not a struct object');
elseif ~all(ismember(p,fieldnames(PROPS)))
    error('wellps:computeParams','Missing property.');    
end

% unpack
aux = num2cell(G.cartDims); 
[I,J,K] = aux{:};


%% Parameters 

PHI  = PROPS.PHI;                     % porosity
KX   = PROPS.KX;                      % permeability x
KY   = PROPS.KY;                      % permeability y
KZ   = PROPS.KZ;                      % permeability z
PHIZ = PHI./(1.0 - PHI);              % normalized porosity

% average permeabilities 
% (It were added several averages here for testing purposes. Note that 
% they are particular cases of the 'generalized mean', or Holder mean)
KA = 1/3*(KX + KY + KY);                        % arithmetic
KG = (KX.*KY.*KZ).^(1/3);                       % geometric
KN = sqrt(KX.^2 + KY.^2 + KZ.^2);               % normalized
KQ = sqrt(1.0/3.0*(KX.^2 + KY.^2 + KZ.^2));     % quadratic
KH = 3.0./(KX.^(-1) + KY.^(-1) + KZ.^(-1));      % harmonic

% reservoir quality index (RQI)
% needs 'isinf' to remove null porosity cells
r = @(var) 0.0314*sqrt(var./PHI);
RQIA = r(KA);  RQIA(isinf(RQIA)) = 0; % arithmetic
RQIG = r(KG);  RQIG(isinf(RQIG)) = 0; % geometric
RQIN = r(KN);  RQIN(isinf(RQIN)) = 0; % normalized
RQIQ = r(KQ);  RQIQ(isinf(RQIQ)) = 0; % quadratic
RQIH = r(KH);  RQIH(isinf(RQIH)) = 0; % harmonic

% flow zone indicator (FZI)
r = @(RQI) RQI./PHIZ;
FZIA = r(RQIA); FZIA(isnan(FZIA)) = 0; % arithmetic
FZIG = r(RQIG); FZIG(isnan(FZIG)) = 0; % geometric
FZIN = r(RQIN); FZIN(isnan(FZIN)) = 0; % normalized
FZIQ = r(RQIQ); FZIQ(isnan(FZIQ)) = 0; % quadratic
FZIH = r(RQIH); FZIH(isnan(FZIH)) = 0; % harmonic

% flow zone indicator * (FZI*)
% Ref.: Paiaman (2015), DOI:10.1002/ente.201500010
FZIAStar = RQIA;
FZIGStar = RQIG;
FZINStar = RQIN;
FZIQStar = RQIQ;
FZIHStar = RQIH;

% discrete rock type (DRT) - neperian log
% Ref. Guo (2005) 
% "Rock Typing as an Effective Tool for Permeability 
% and Water-Saturation Modeling: A Case Study in a Clastic Reservoir 
% in the Oriente Basin"
r = @(FZI) round(2*log(FZI) + 10.6);
DRTA_LN = r(FZIA); DRTA_LN(isinf(DRTA_LN)) = 0; % arithmetic
DRTG_LN = r(FZIG); DRTG_LN(isinf(DRTG_LN)) = 0; % geometric
DRTN_LN = r(FZIN); DRTN_LN(isinf(DRTN_LN)) = 0; % normalized
DRTQ_LN = r(FZIQ); DRTQ_LN(isinf(DRTQ_LN)) = 0; % quadratic
DRTH_LN = r(FZIH); DRTH_LN(isinf(DRTH_LN)) = 0; % harmonic

% DRT*, for FZI*
DRTAStar_LN = r(FZIAStar); DRTAStar_LN(isinf(DRTAStar_LN)) = 0; % arithmetic
DRTGStar_LN = r(FZIGStar); DRTGStar_LN(isinf(DRTGStar_LN)) = 0; % geometric
DRTNStar_LN = r(FZINStar); DRTNStar_LN(isinf(DRTNStar_LN)) = 0; % normalized
DRTQStar_LN = r(FZIQStar); DRTQStar_LN(isinf(DRTQStar_LN)) = 0; % quadratic
DRTHStar_LN = r(FZIHStar); DRTHStar_LN(isinf(DRTHStar_LN)) = 0; % harmonic

% discrete rock type (DRT) - base-10 log
r = @(FZI) round(2*log10(FZI) + 10.6);
DRTA_LOG10 = r(FZIA); DRTA_LOG10(isinf(DRTA_LOG10)) = 0; % arithmetic
DRTG_LOG10 = r(FZIG); DRTG_LOG10(isinf(DRTG_LOG10)) = 0; % geometric
DRTN_LOG10 = r(FZIN); DRTN_LOG10(isinf(DRTN_LOG10)) = 0; % normalized
DRTQ_LOG10 = r(FZIQ); DRTQ_LOG10(isinf(DRTQ_LOG10)) = 0; % quadratic
DRTH_LOG10 = r(FZIH); DRTH_LOG10(isinf(DRTH_LOG10)) = 0; % harmonic

% DRT*, for FZI*
DRTAStar_LOG10 = r(FZIAStar); DRTAStar_LOG10(isinf(DRTAStar_LOG10)) = 0; % arithmetic
DRTGStar_LOG10 = r(FZIGStar); DRTGStar_LOG10(isinf(DRTGStar_LOG10)) = 0; % geometric
DRTNStar_LOG10 = r(FZINStar); DRTNStar_LOG10(isinf(DRTNStar_LOG10)) = 0; % normalized
DRTQStar_LOG10 = r(FZIQStar); DRTQStar_LOG10(isinf(DRTQStar_LOG10)) = 0; % quadratic
DRTHStar_LOG10 = r(FZIHStar); DRTHStar_LOG10(isinf(DRTHStar_LOG10)) = 0; % harmonic

%{ 
    REMARK
    ======

    Please, note that the arrays that compute Log values might produce -Inf
    /+Inf indeterminacies because there are points with null porosity. 
    This was observed afterwards and a forced conversion to zero was
    inserted here.
    
    User must be aware of this, since the outcome arrays will not bring
    arrays with the original values. For the purposes of this project, null
    porosity cells have no sense, as with keeping arrays with -/+ Inf values. 
    The 0 entries, however, must be therein to keep the coherent dimensions
    for final plotting and addressing. 

    Furthermore, we have not chosen sparse allocation because
    incompatibilities with MRST might arise, since the functions would have
    to receive sparse objects as input arguments. 

    Below, we make the 0 filtering in the same manner we did for DRT arrays. 

%}

% base-10 logs
Log10PHIZ = log10(PHIZ);  Log10PHIZ(isinf(Log10PHIZ)) = 0;   
Log10RQIA = log10(RQIA);  Log10RQIA(isinf(Log10RQIA)) = 0;   % arithmetic                 
Log10RQIG = log10(RQIG);  Log10RQIG(isinf(Log10RQIG)) = 0;   % geometric
Log10RQIN = log10(RQIN);  Log10RQIN(isinf(Log10RQIN)) = 0;   % normalized
Log10RQIQ = log10(RQIQ);  Log10RQIQ(isinf(Log10RQIQ)) = 0;   % quadratic
Log10RQIH = log10(RQIH);  Log10RQIH(isinf(Log10RQIH)) = 0;   % harmonic

% base-e logs
LNPHIZ = log(PHIZ);  LNPHIZ(isinf(LNPHIZ)) = 0;    
LNRQIA = log(RQIA);  LNRQIA(isinf(LNRQIA)) = 0;  % arithmetic                 
LNRQIG = log(RQIG);  LNRQIG(isinf(LNRQIG)) = 0;  % geometric
LNRQIN = log(RQIN);  LNRQIN(isinf(LNRQIN)) = 0;  % normalized
LNRQIQ = log(RQIQ);  LNRQIQ(isinf(LNRQIQ)) = 0;  % quadratic
LNRQIH = log(RQIH);  LNRQIH(isinf(LNRQIH)) = 0;  % harmonic

% for FZI*, we need log(0.0314*sqrt(KN)) and log(sqrt(PHI))
% # base-10 logs
Log10FZIStar_SQRTPHI = log10(sqrt(PHI)); 
Log10FZIStarA_SQRTK = log10(0.0314*sqrt(KA)); % arithmetic                 
Log10FZIStarG_SQRTK = log10(0.0314*sqrt(KG)); % geometric
Log10FZIStarN_SQRTK = log10(0.0314*sqrt(KN)); % normalized
Log10FZIStarQ_SQRTK = log10(0.0314*sqrt(KQ)); % quadratic
Log10FZIStarH_SQRTK = log10(0.0314*sqrt(KH)); % hamonic

Log10FZIStar_SQRTPHI(isinf(Log10FZIStar_SQRTPHI)) = 0;
Log10FZIStarA_SQRTK(isinf(Log10FZIStarA_SQRTK)) = 0;
Log10FZIStarG_SQRTK(isinf(Log10FZIStarG_SQRTK)) = 0;
Log10FZIStarN_SQRTK(isinf(Log10FZIStarN_SQRTK)) = 0;
Log10FZIStarQ_SQRTK(isinf(Log10FZIStarQ_SQRTK)) = 0;
Log10FZIStarH_SQRTK(isinf(Log10FZIStarH_SQRTK)) = 0;


% # base-e logs
LNFZIStar_SQRTPHI = log(sqrt(PHI)); 
LNFZIStarA_SQRTK = log(0.0314*sqrt(KA)); % arithmetic                 
LNFZIStarG_SQRTK = log(0.0314*sqrt(KG)); % geometric
LNFZIStarN_SQRTK = log(0.0314*sqrt(KN)); % normalized
LNFZIStarQ_SQRTK = log(0.0314*sqrt(KQ)); % quadratic
LNFZIStarH_SQRTK = log(0.0314*sqrt(KH)); % hamonic

LNFZIStar_SQRTPHI(isinf(LNFZIStar_SQRTPHI)) = 0;
LNFZIStarA_SQRTK(isinf(LNFZIStarA_SQRTK)) = 0;
LNFZIStarG_SQRTK(isinf(LNFZIStarG_SQRTK)) = 0;
LNFZIStarN_SQRTK(isinf(LNFZIStarN_SQRTK)) = 0;
LNFZIStarQ_SQRTK(isinf(LNFZIStarQ_SQRTK)) = 0;
LNFZIStarH_SQRTK(isinf(LNFZIStarH_SQRTK)) = 0;


% all params
p = {PHI, KX, KY, KZ,  PHIZ, KA, KG, KN, KQ, KH,                      ...
     RQIA, RQIG, RQIN, RQIQ, RQIH, FZIA, FZIG, FZIN, FZIQ, FZIH,      ...
     FZIAStar, FZIGStar, FZINStar, FZIQStar, FZIHStar,                ...
     DRTA_LN, DRTG_LN, DRTN_LN, DRTQ_LN, DRTH_LN,                     ...
     DRTA_LOG10, DRTG_LOG10, DRTN_LOG10, DRTQ_LOG10, DRTH_LOG10,      ...
     DRTAStar_LN,                                                     ...
     DRTGStar_LN,                                                     ...
     DRTNStar_LN,                                                     ...
     DRTQStar_LN,                                                     ...
     DRTHStar_LN,                                                     ...
     DRTAStar_LOG10,                                                  ...
     DRTGStar_LOG10,                                                  ...
     DRTNStar_LOG10,                                                  ...
     DRTQStar_LOG10,                                                  ...
     DRTHStar_LOG10,                                                  ...
     Log10PHIZ, Log10RQIA, Log10RQIG, Log10RQIN, Log10RQIQ, Log10RQIH,...
     LNPHIZ, LNRQIA, LNRQIG, LNRQIN, LNRQIQ, LNRQIH,                  ...
     Log10FZIStarA_SQRTK,                                             ...
     Log10FZIStarG_SQRTK,                                             ...
     Log10FZIStarN_SQRTK,                                             ...
     Log10FZIStarQ_SQRTK,                                             ...
     Log10FZIStarH_SQRTK,                                             ...
     Log10FZIStar_SQRTPHI,                                            ...
     LNFZIStarA_SQRTK,                                                ...
     LNFZIStarG_SQRTK,                                                ...
     LNFZIStarN_SQRTK,                                                ...
     LNFZIStarQ_SQRTK,                                                ...
     LNFZIStarH_SQRTK,                                                ...
     LNFZIStar_SQRTPHI};
     
   
% field names
s = {'PHI', 'KX', 'KY', 'KZ', 'PHIZ', 'KA', 'KG', 'KN', 'KQ', 'KH',   ...
     'RQIA', 'RQIG', 'RQIN', 'RQIQ', 'RQIH',                          ...
     'FZIA', 'FZIG', 'FZIN', 'FZIQ', 'FZIH',                          ...
     'FZIAStar', 'FZIGStar', 'FZINStar', 'FZIQStar', 'FZIHStar',      ...
     'DRTA_LN', 'DRTG_LN', 'DRTN_LN', 'DRTQ_LN', 'DRTH_LN',           ...
     'DRTA_LOG10',                                                    ...
     'DRTG_LOG10',                                                    ...
     'DRTN_LOG10',                                                    ...
     'DRTQ_LOG10',                                                    ...
     'DRTH_LOG10',                                                    ...
     'DRTAStar_LN',                                                   ...
     'DRTGStar_LN',                                                   ...
     'DRTNStar_LN',                                                   ...
     'DRTQStar_LN',                                                   ...
     'DRTHStar_LN',                                                   ...
     'DRTAStar_LOG10',                                                ...
     'DRTGStar_LOG10',                                                ...
     'DRTNStar_LOG10',                                                ...
     'DRTQStar_LOG10',                                                ...
     'DRTHStar_LOG10',                                                ...
     'Log10PHIZ',                                                     ... 
     'Log10RQIA',                                                     ... 
     'Log10RQIG',                                                     ... 
     'Log10RQIN',                                                     ... 
     'Log10RQIQ',                                                     ... 
     'Log10RQIH',                                                     ...    
     'LNPHIZ', 'LNRQIA', 'LNRQIG', 'LNRQIN', 'LNRQIQ', 'LNRQIH',      ...
     'Log10FZIStarA_SQRTK',                                           ...
     'Log10FZIStarG_SQRTK',                                           ...
     'Log10FZIStarN_SQRTK',                                           ...
     'Log10FZIStarQ_SQRTK',                                           ...
     'Log10FZIStarH_SQRTK',                                           ...
     'Log10FZIStar_SQRTPHI',                                          ...
     'LNFZIStarA_SQRTK',                                              ...
     'LNFZIStarG_SQRTK',                                              ...
     'LNFZIStarN_SQRTK',                                              ...
     'LNFZIStarQ_SQRTK',                                              ...
     'LNFZIStarH_SQRTK',                                              ...
     'LNFZIStar_SQRTPHI'};

% reshape all
fp = @(param) reshape(param,[I,J,K]);
for i = 1:numel(p), p{i} = fp(p{i});  end

% unfold into struct
P = cell2struct(p,s,2);
