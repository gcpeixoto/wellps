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

%plotCellData(G,PROPS.SO(PROPS.ACTNUM == 1))

%% Mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;


%% Oil-saturated cells 
%idso = find(PROPS.SO(PROPS.ACTNUM == 1) > 0);
%id = find(Ind(idso));
%subG = extractSubgrid(G,id);
%plotCellData(subG,find(Ind(PROPS.SO(id))));

% parameters
%getPso = @(x) x(idso); 
%Pso = structfun(getPso,P,'UniformOutput',false);
%Sso = getStats(d,Pso,{'DRTH_LN'},'n');



