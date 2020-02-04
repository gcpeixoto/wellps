%% TUTORIAL: main script that demonstrates the WELLPS walkthrough

mrstVerbose off  % turn on verbose

%% Mounting 

% class instantiation 
d = DirManager(); 

%% Grid reading
f = fullfile(d.getBenchMarksDir,'unisim-I-D','eclipse','UNISIM_I_D_ECLIPSE.DATA');

[G,PROPS] = buildModel(f);

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
% getStats is useful to show a short summary of frequencies and values 
% of the variables computed in P

% field statistics for chosen properties
S = getStats(d,P,{'DRTN_LN','DRTN_LOG10'},'n');

%% Find DRT connections
% Compute DRT connections based on list chosen from observation of
% statistical information in S. It is recommended to expurge DRT = 0.
% This will produce a big structure containing several structures, each per
% DRT value. The flag 'tocsv' allows us also to save this information in 
% .csv files.

% number of significant cells
nofs = 30;

drtlist = S{1}(5:6,1);
drtSt = findDRTConnections(d,drtlist, P, 'geometric','ln',nofs,'y', 1);

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
opt.nofs = nofs;    
% linear regression slope +/- tolerance.
opt.seps = 0.05;
% minimum R2 coefficient tolerance.
opt.R2min = 0.9;

% compute
[metricsSt,linregrSt] = computeDRTGraphMetrics(opt,drtSt);


%% Metrics analytics 
% This method will export lots of .csv files ready to be handled in terms
% of data analytics relating to graph centrality metrics computed for each 
% cluster.  
metricsAnalyzer([13,14],[1,2],'geometric','ln');

%% Process cluster fit 
% This method will compute the fit parameters for given DRT,
% cluster(s) and rotation angle(s) based on a best-fit ellipsoid
% for the cluster shape. These parameters are useful to build 
% 5-spot nonuniform well patterns. See the functioin documentation to 
% understand what is going on behind the scenes. 
clusterFitSt = processClusterFit(G,[13,14],1:4,'geometric','ln',[0,pi/2]);

%% Build a nonuniform 5-spot well pattern for a given cluster
% Here, we form the nonuniform standard 5-spot as follows:
%
% - producer well: column of all (vertical) cells neighbour to the maximum
%                  closeness centrality cell 
% - injection wells: column of all (vertical) cells neighbour to each of 
%                    the 4 cluster-fit cells +X,-X,+Y,-Y for the cluster
%
% REMARK: the z-range considered is that one of the cluster.

% injectors 
% +X
col_inj_1 = clusterFitSt.DRT13.C1.Pattern1.colNeighsX1;
                  
% -X
col_inj_2 = clusterFitSt.DRT13.C1.Pattern1.colNeighsX2;
         
% +Y
col_inj_3 = clusterFitSt.DRT13.C1.Pattern1.colNeighsY1;
                  
% -Y
col_inj_4 = clusterFitSt.DRT13.C1.Pattern1.colNeighsY2;

% producer 
col_prod = clusterFitSt.DRT13.C1.Pattern1.colNeighsMaxC;

% getting cell indices
icol_prod = sub2ind(G.cartDims,col_prod(:,1),col_prod(:,2),col_prod(:,3));
icol_prod = Ind(icol_prod);

icol_inj_1 = sub2ind(G.cartDims,col_inj_1(:,1),col_inj_1(:,2),col_inj_1(:,3));
icol_inj_1 = Ind(icol_inj_1);

icol_inj_2 = sub2ind(G.cartDims,col_inj_2(:,1),col_inj_2(:,2),col_inj_2(:,3));
icol_inj_2 = Ind(icol_inj_2);

icol_inj_3 = sub2ind(G.cartDims,col_inj_3(:,1),col_inj_3(:,2),col_inj_3(:,3));
icol_inj_3 = Ind(icol_inj_3);

icol_inj_4 = sub2ind(G.cartDims,col_inj_4(:,1),col_inj_4(:,2),col_inj_4(:,3));
icol_inj_4 = Ind(icol_inj_4);

% plot 5-spot columns 
figure 
plotGrid(G, Ind(drtSt.DRT13.compVoxelInds{1}),'FaceColor',[0.6,0.6,0.6], ...
    'FaceAlpha',0.05, 'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.1)

% producer
plotGrid(G,icol_prod(~isnan(icol_prod)),'FaceColor','r','EdgeColor','k')

% injector 1
plotGrid(G,icol_inj_1(~isnan(icol_inj_1)),'FaceColor','b','EdgeColor','k')

% injector 2
plotGrid(G,icol_inj_2(~isnan(icol_inj_2)),'FaceColor','b','EdgeColor','k')

% injector 3
plotGrid(G,icol_inj_3(~isnan(icol_inj_3)),'FaceColor','b','EdgeColor','k')

% injector 4
plotGrid(G,icol_inj_4(~isnan(icol_inj_4)),'FaceColor','b','EdgeColor','k')
axis off vis3d
