function [J,G,PROPS,active] = computeProdProxy(f,proxy,bhp)
%COMPUTEPRODPROXY Compute productivity potential functions over 
%                 a given grid.
%
%   PARAMETERS:
%
%       f       -  path to ECLIPSE deck file (char)
%       proxy   -  proxy function to compute (char). 
%                  Currently implemented methods are:
%                   'rqi', 'rqip', 'kharghoria', 'liu', and 'rqipd'. 
%                  See explanation below.
%
%   RETURNS:
%
%       J       - productivity potential function (3D array)
%       G       - processed Eclipse grid (struct)
%       PROPS   - original Eclipse grid (struct)
%       active  - active cell indices (1D array)
%
%
%
% REMARK: 'liu' and 'rqipd' MUST BE TESTED FURTHER.

%% Grid reading and processing
% We are unable to call the function 'lib/buildModel' because we need
% nonstandard fields forcibly inserted into MRST:readGRDECL function.
% Then we make adaptions to retrieve oil saturation and initial reservoir 
% pressure. Read the REMARK in the case study main M-file.

% get grid
G = readGRDECL(f); 

% checking
assert(isfield(G,'SO') && isfield(G,'PRESSURE'),...
    ['Nonstandards fields ''SO'' and ''PRESSURE'' ', ... 
     'are missing in this grid. Have you read the initial REMARK about this?'])    

% create an auxiliary 'PROPS' struct to reuse G
PROPS.PHI = G.PORO;
PROPS.KX  = G.PERMX;
PROPS.KY  = G.PERMY;
PROPS.KZ  = G.PERMZ;
PROPS.ACTNUM  = G.ACTNUM;
PROPS.SO  = G.SO; % oil sat
PROPS.PRESS0  = G.PRESSURE; % init press

% grid processing
G = processGRDECL(G);
G = computeGeometry(G);

% compute field attributes (poro, perm, RQI, etc.)
P = computeParams(G,PROPS); 

% active cells
active = find(PROPS.ACTNUM == true);

%% Productivity potential function
% By following the literature, we denote the PP function by 'J'.
% Here, the J-function is taken cell by cell and stands for local 
% productivity.
%
% We model J as a space-time function J(c,t) where c is the cell location
% and t is the time instant. However, in this script, we consider t = 0.
%
% Below, we compute J through different methods as follows:
%
%   - 'rqi':  defines J = RQIn*SO = J1, 
%
%      where RQIn is the normalized RQI and SO is the oil saturation.
%      Under development by Gustavo @LaMEP.
%
%
%   - 'rqip': defines J2 = J1*Pn, 
%
%      where Pn is the normalized pressure (P - Pmin)/(Pmax - Pmin). 
%      Under development by Gustavo @LaMEP.
%
%   - 'kharghoria': defines J3 = K*SO*PHI, 
%
%      where K is the absolute permeability and PHI is the porosity.
%      Here, we assume K = KNn so that J3 is ALSO NORMALIZED.
%      Originally introduced by (Kharghoria et al., 2003)                          
%
%   - 'liu': defines J4 = (SO - SOres)*(Po - Pomin)*log(K)*log(RDIST(c,B))
%
%      where SOres is the residual oil saturation, Po is oil-phase pressure
%      K is the absolute permeability and RDIST is the distance
%      between the grid cell c and the nearest boundary.
%      Here, we consider SOres = 0, Po = usual hydrostatic pressure and
%      B is formed by all cells that have an exterior face. 
%      Besides, J4 is ALSO NORMALIZED.
%      Originally introduced by (Liu et al., 2006)  
%
%   - 'rqipd': defines J5 = J2*log(RDIST(c,B))
%
%      This adds the distance to the boundary to our approach RQI*SO*P.

% -- Required variables

% norm-based RQI = sqrt(KN/PHI)
RQIn = P.RQIN;

%  norm-based permeability
KNn = P.KN;

% pressure differential 
% (resevoir pressure - bottom hole pressure)
% Here, we assume that all wells in the field operate at the same fixed
% BHP.
PRESSn = PROPS.PRESS0 - bhp;
     
% reshape is necessary to operate connections over all the grid, 
% not only over active cells.
PHI_r = reshape(PROPS.PHI,G.cartDims);
PRESSn_r = reshape(PRESSn,G.cartDims);    
SO_r = reshape(PROPS.SO,G.cartDims);

% -- J calculation
switch proxy
    
    case 'rqi'            
        J = RQIn.*SO_r;            
    
    case 'rqip'        
        J = RQIn.*SO_r.*PRESSn_r;                                
            
    case 'kharghoria'                
        J = PHI_r.*SO_r.*KNn;        
                    
    case 'liu'                
        Rn = computeDistBoundary(G);                
        J = SO_r.*PRESSn_r.*log(KN).*log(Rn);
        
    case 'rqipd'
        Rn = computeDistBoundary(G);
        J = RQIn.*SO_r.*PRESSn_r.*log(Rn);      
        
    otherwise
        error('Productivity proxy function not implemented!')
                    
end

% set all NANs and INFs entries to zero potential.
J(isnan(J)) = 0;
J(isinf(J)) = 0;

% normalization
J = J./max(J);

end