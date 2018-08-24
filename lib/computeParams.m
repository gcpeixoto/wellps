function P = computeParams(G,PROPS)
% Compute petrophysical and analytical parameters required to
% reservoir modelling

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

% discrete rock type (DRT)
r = @(FZI) round(2*log(FZI) + 10.6);
DRTA = r(FZIA); DRTA(isinf(DRTA)) = 0; % arithmetic
DRTG = r(FZIG); DRTG(isinf(DRTG)) = 0; % geometric
DRTN = r(FZIN); DRTN(isinf(DRTN)) = 0; % normalized

% base-10 logs
LogPHIZ = log10(PHIZ);    
LogRQIA = log10(RQIA); % arithmetic                 
LogRQIG = log10(RQIG); % geometric
LogRQIN = log10(RQIN); % normalized



% all params
p = {PHI, KX, KY, KZ , PHIZ,    ...
     KA, KG, KN,                ...
     RQIA, RQIG, RQIN,          ...
     FZIA, FZIG, FZIN,          ...
     DRTA, DRTG, DRTN,          ...
     LogPHIZ,                   ...
     LogRQIA, LogRQIG, LogRQIN};
   
% field names
s = {'PHI','KX','KY','KZ','PHIZ',   ...
     'KA','KG','KN',                ...
     'RQIA','RQIG','RQIN',          ...
     'FZIA','FZIG','FZIN',          ...
     'DRTA','DRTG','DRTN',          ...
     'LogPHIZ',                     ...
     'LogRQIA','LogRQIG','LogRQIN'};

% reshape all
fp = @(param) reshape(param,[I,J,K]);
for i = 1:numel(p), p{i} = fp(p{i});  end

% unfold into struct
P = cell2struct(p,s,2);
