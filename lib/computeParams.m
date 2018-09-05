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
p = {'PHI','KX','KY','KZ','ACTNUM'};
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
KA = 1/3*(KX + KY + KY);              % arithmetic
KG = (KX.*KY.*KZ).^(1/3);             % geometric
KN = sqrt(KX.^2 + KY.^2 + KZ.^2);     % normalized

% reservoir quality index (RQI)
% needs 'isinf' to remove null porosity cells
r = @(var) 0.0314*sqrt(var./PHI);
RQIA = r(KA);  RQIA(isinf(RQIA)) = 0; % arithmetic
RQIG = r(KG);  RQIG(isinf(RQIG)) = 0; % geometric
RQIN = r(KN);  RQIN(isinf(RQIN)) = 0; % normalized

% flow zone indicator (FZI)
r = @(RQI) RQI./PHIZ;
FZIA = r(RQIA); FZIA(isnan(FZIA)) = 0; % arithmetic
FZIG = r(RQIG); FZIG(isnan(FZIG)) = 0; % geometric
FZIN = r(RQIN); FZIN(isnan(FZIN)) = 0; % normalized

% discrete rock type (DRT) - neperian log
r = @(FZI) round(2*log(FZI) + 10.6);
DRTA_LN = r(FZIA); DRTA_LN(isinf(DRTA_LN)) = 0; % arithmetic
DRTG_LN = r(FZIG); DRTG_LN(isinf(DRTG_LN)) = 0; % geometric
DRTN_LN = r(FZIN); DRTN_LN(isinf(DRTN_LN)) = 0; % normalized

% discrete rock type (DRT) - neperian log
r = @(FZI) round(2*log10(FZI) + 10.6);
DRTA_LOG10 = r(FZIA); DRTA_LOG10(isinf(DRTA_LOG10)) = 0; % arithmetic
DRTG_LOG10 = r(FZIG); DRTG_LOG10(isinf(DRTG_LOG10)) = 0; % geometric
DRTN_LOG10 = r(FZIN); DRTN_LOG10(isinf(DRTN_LOG10)) = 0; % normalized

% base-10 logs
Log10PHIZ = log10(PHIZ);    
Log10RQIA = log10(RQIA); % arithmetic                 
Log10RQIG = log10(RQIG); % geometric
Log10RQIN = log10(RQIN); % normalized

% base-2 logs
LNPHIZ = log(PHIZ);    
LNRQIA = log(RQIA); % arithmetic                 
LNRQIG = log(RQIG); % geometric
LNRQIN = log(RQIN); % normalized


% all params
p = {PHI, KX, KY, KZ , PHIZ,                     ...
     KA, KG, KN,                                 ...
     RQIA, RQIG, RQIN,                           ...
     FZIA, FZIG, FZIN,                           ...
     DRTA_LN, DRTG_LN, DRTN_LN,                  ...
     DRTA_LOG10, DRTG_LOG10, DRTN_LOG10,         ...
     Log10PHIZ, Log10RQIA, Log10RQIG, Log10RQIN, ...
     LNPHIZ, LNRQIA, LNRQIG, LNRQIN};
   
% field names
s = {'PHI','KX','KY','KZ','PHIZ',                      ...
     'KA','KG','KN',                                   ...
     'RQIA','RQIG','RQIN',                             ...
     'FZIA','FZIG','FZIN',                             ...
     'DRTA_LN', 'DRTG_LN', 'DRTN_LN',                  ...
     'DRTA_LOG10', 'DRTG_LOG10', 'DRTN_LOG10',         ...
     'Log10PHIZ', 'Log10RQIA','Log10RQIG','Log10RQIN', ...
     'LNPHIZ', 'LNRQIA','LNRQIG','LNRQIN'};

% reshape all
fp = @(param) reshape(param,[I,J,K]);
for i = 1:numel(p), p{i} = fp(p{i});  end

% unfold into struct
P = cell2struct(p,s,2);
