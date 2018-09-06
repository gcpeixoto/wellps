%% TUTORIAL 

mrstVerbose off  % turn on verbose

%% Mounting 

% class instantiations 
d = DirManager(); 

d.mountDir();   % mounts standard directory tree

%% Grid reading
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');

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
drtSt = findDRTConnections(drtlist, P, 'geometric','ln',500,'y', 1);

%% Plot clusters from DRT-connected cells
% To plot DRT-connected clusters, we need to take two steps: i) to use the
% list of cell indices relative to the cluster you want to plot; ii) to do
% an 'inverse mapping'. 
% For i), it suffices to collect one among the lists contained inside 
% the 'drtSt' structure; for ii), we need to follow a few steps
%

% This makes the 'inverse mapping', in the sense of filling the positions of 
% Ind with the original indices of the grid before being processed by 
% 'processGRDECL'. The array 'indexMap' has the original indices to which 
% G.cells are identified.
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

% plots the full UNISIM in dimmed grey 
figure
plotGrid(G, 1:G.cells.num,'FaceColor',[0.6,0.6,0.6], ...
    'FaceAlpha',0.05, 'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.)

% plots clusters 5 and 8 for DRT = 13.
% We assume here DRT13 was included in 'drtlist' above and computed! 
plotGrid(G, Ind(drtSt.DRT13.compVoxelInds{5}),... 
    'FaceColor','c','EdgeColor','k')
plotGrid(G, Ind(drtSt.DRT13.compVoxelInds{8}),...
    'FaceColor','r','EdgeColor','k')
axis off vis3d

%% Computing graph metrics 
% Here, we choose parameters to compute the graph metrics over all
% the clusters previously computed with findDRTConnections 

% number of significant cells to be considered. Clusters with less 'nofs'
% cells are ignored
opt.nofs = 100;    
% linear regression slope +/- tolerance.
opt.seps = 0.05;
% minimum R2 coefficient tolerance.
opt.R2min = 0.9;

% compute
[metricsSt,linregrSt] = computeDRTGraphMetrics(opt,drtSt);

