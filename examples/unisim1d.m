%% TUTORIAL 

mrstVerbose on  % turn on verbose

%% Mounting 

% class instantiations 
d = DirManager(); 

d.mountDir();   % mounts standard directory tree

%% Grid reading
%[G,GAll,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');
[G,GAll,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING_NOPINCH.DATA');

%% Plot properties 
% 'buildModel' delivers a processed grid structure G.
% In ECLIPSE decks, MRST bypasses cells whose ACTNUM property
% is null. Hence, the number of cells in G.cells
% causes porosity and other properties in PROPS to have 
% higher dimension. To plot properties, we need to ignore
% ACTNUM. 
on = PROPS.ACTNUM == 1;

% plot porosity
%plotCellData(G,PROPS.PHI(on));

% plot permeability x
%figure
%plotCellData(G,PROPS.KX(on));

%% Compute required parameters
% Because of ACTNUM, we need to work with all the cells, thus 
% computing the required parameters even in the points where ACTNUM = 0. 
% Note that the total number of cells is the product of the 
% princiapal discretization 
ncells = prod(G.cartDims) == numel(PROPS.ACTNUM); % returns 'true'

% structure of several parameters (RQI, FZI, PHIZ, etc.)
P = computeParams(G,PROPS);

%% Statistics for field variables
% printStats is useful to show a short summary of frequencies and values 
% of the variables computed in P

% field statistics for chosen properties
S = printStats(P,{'DRTN_LN','DRTN_LOG10'},'n');

%% Find DRT connections
% Compute DRT connections based on list chosen from observation of
% statistical information in S. It is recommended to expurge DRT = 0.
% This will produce a big structure containing several structures, each per
% DRT value. The flag 'tocsv' allows us also to save this information in 
% .csv files.
drtlist = S{1}(5:6,1);
drtSt = findDRTConnections(drtlist, P.DRTN_LN,500,'n', 1);

