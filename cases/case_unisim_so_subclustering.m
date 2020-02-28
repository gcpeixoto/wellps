%% HFU sub-clustering 
% Based on oil-saturated cells
%
% To run this case, the MRST's function 'readGRDECL' was modified to
% include the keyword 'SO', which is a user-defined keyword for oil
% saturation.
%
% ---> MODIFICATION TO BE MADE IN 'readGRDECL' (beetween ** .. **)
%
%           case {'PORO',                         ...
%                'PERMX' , 'PERMXY', 'PERMXZ',   ...
%                'PERMYX', 'PERMY' , 'PERMYZ',   ...
%                'PERMZX', 'PERMZY', 'PERMZ' ,   ...
%                'PERMH',                        ...
%                'ACTNUM', 'SATNUM', 'ROCKTYPE', ...
%                'MULTX' , 'MULTX-',             ...
%                'MULTY' , 'MULTY-',             ...
%                'MULTZ' , 'MULTZ-',             ...
%                'NTG'   , 'VSH'   , **'SO'**    ...
%                }
%
% Moreover, for this case, we used a duplicata file 'UNISIM_....SO.DATA'
% which has a INCLUDE call to the file 'SO.INC' for the distribution 
% of oil saturation obtained by the standard relative permeability curves 
% of the model. 

case_name = 'unisim1-subclustering';  

%% Loading grid 

d = DirManager(); 

f = fullfile(d.getBenchMarksDir,'unisim-I-D','eclipse','UNISIM_I_D_ECLIPSE_SO.DATA');

% Unable to call buildModel because SO is a nonstandard field
G = readGRDECL(f);
PROPS.PHI = G.PORO;
PROPS.KX  = G.PERMX;
PROPS.KY  = G.PERMY;
PROPS.KZ  = G.PERMZ;
PROPS.ACTNUM  = G.ACTNUM;
PROPS.SO  = G.SO;
G = processGRDECL(G);

%% Compute parameters 
P = computeParams(G,PROPS); 

%plotCellData(G,PROPS.SO(on))

%% Mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;


%% Oil-saturated cells 
idso = find(PROPS.SO > 0);
SOP = PROPS.SO(idso);

% zero oil sat
%idsoz = find(PROPS.SO == 0);
%SOZ = PROPS.SO(idsoz);

idsoloc = Ind(idso);

subG = extractSubgrid(G,idsoloc);
plotCellData(subG,SOP,'EdgeColor','None')

% parameters
getPso = @(x) x(idso); 
Pso = structfun(getPso,P,'UniformOutput',false);
Sso = getStats(d,Pso,{'DRTH_LN'},'n');


%% DRT choice 

% get all DRTs 
drtlist = Sso{1}(:,1); 
drtlist = drtlist(drtlist > 0);

% compute HFUs by DRT
nofsc = 10;
drtSt = findDRTConnections(d,drtlist, Pso, 'harmonic','ln',nofsc,'n', 1);



