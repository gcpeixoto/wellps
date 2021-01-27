%% EXAMPLE: centrality metrics from NETWORKX
%
% This example shows that betweeness is incorrectly computed by SNAP, since 
% negative values are found. Closeness, on the other hand, seems correct.
% Degree centrality cannot be directly compared because in NETWORKX it is
% fractional.
% 
% That is why we will use NETWORKX as standard to compute betweeness.

mrstVerbose off  % turn on verbose

case_name = 'ex4';

%% Mounting 

% class instantiation 
d = DirManager(); 

%% Loading grid 
f = fullfile(d.getBenchMarksDir,'psy','eclipse','PSY.grdecl');

%% Productivity Potential Index 
% Compute PPI through one of the available methods:
m = {'rqi', 'rqip', 'kharghoria'};
[J,G,PROPS,active] = computeProdProxy(f,m{2});

%% Productivity unit classes
% Computation of the PUCs (productivity unit class).
% For a given reservoir cell w and proxy function J, J : w -> [0,1]. Then
% we define
%
%          { 1, if 0.0 < F(w) <  d1
% PUC(w) = { 2, if d1 <= F(w) <  d2
%          { 3, if d2 <= F(w) <  d3
%          { ...
%          { n, if dn <= F(w) <= max(F)
%
%
% Therewith, we define a PU of class X as each connected cluster 
% whose all cells have PUC(w) = X.

% compute discrete function PUC
[PUC,nclasses,delta,div] = computePUC(J,active,'3');

% list of all PUC values except zero
puc = 1:nclasses;

% 6-neighbor connectivity criterion to find clusters from the PUC values
nofsc = 50; % only for .csv
pucSt = findConnectionsByPUC(d,1:nclasses,PUC,nofsc,'y',1);

%% Computing graph metrics 
% Here, we choose parameters to compute the graph metrics over all
% the clusters previously computed with findDRTConnections 

% number of significant cells to be considered. Clusters with less 'nofsc'
% cells are ignored
opt.nofsc = nofsc;    

% Compute with NETWORKX
Mf = computePUCGraphMetricsNetworkx(opt,pucSt);

%% Metrics analysis 

% options
analytics.loaddir = d.getMatDir;
analytics.puclist = puc;
analytics.savedir = fullfile(d.getCsvDir,case_name);
analytics.filefreq = true;
analytics.fileminmax = true;
analytics.filemetrics = true;

% post-processed files
dataDir = metricsAnalysisPUC(analytics);

%% Plotting

% getting cluster ID to plot
% The maximum value is 'pucSt.(mf).allNComps'

nc = 2; % <--- choose here!!!

mf = fieldnames(Mf); 
fprintf('PUC values available: %s\n',char(mf{:}))

mf = mf{1}; % <--- choose here!!!
fprintf('Choosing: %s\n',mf)

mf_max = pucSt.(mf).allNComps;
if nc > mf_max    
    nc = mf_max;
    fprintf('Cluster ID set to %d',nc);
end

%% Inverse map to recover current logical indices at the original (unprocessed) grid
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

%% Plot the background grid 
% This will plot the cluster region in dimmed gray.
% You may control the colors and visual.

figure
plotGrid(G, Ind(pucSt.(mf).compVoxelInds{nc}),...
    'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.05, ...
    'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.5);
axis off, view([82,8])
clf

%% Plot the centroids 
% To plot the cell centroids in space we need to get 
% the global indices of each cluster cell according to 'indexMap'.
% Then, we get the global stored indices from 'compVoxelInds'
% and store them as global indices in the new map.
globInd = Ind(pucSt.(mf).compVoxelInds{nc});


% centroid coordinates of cluster cells
xc = G.cells.centroids(globInd,1);
yc = G.cells.centroids(globInd,2);
zc = G.cells.centroids(globInd,3);
scatter3(xc,yc,zc,30,[0.6,0.2,0.1],'filled')

%% Plot the 3D graph
% First of all, we need to load the adjacency matrix of the cluster, which
% gives us the (local) indices of connected cells

% load M for the chosen
load(Mf.(mf)); 

adj = M.adjMatrix{nc};
[Iloc,Jloc] = find(adj==1);

% Once the indices are found, we create arrays to prepare the 'edge'
% plotting by choosing the source and target cells.

% source centroids
xcsrc = G.cells.centroids(globInd(Iloc),1); 
ycsrc = G.cells.centroids(globInd(Iloc),2); 
zcsrc = G.cells.centroids(globInd(Iloc),3); 

% target centroids
xctgt = G.cells.centroids(globInd(Jloc),1); 
yctgt = G.cells.centroids(globInd(Jloc),2); 
zctgt = G.cells.centroids(globInd(Jloc),3); 

hold on
% This will plot connected lines to mimic the graph edges
for i = 1:numel(Iloc)
    lin = line([xcsrc(i),xctgt(i)], ... 
               [ycsrc(i),yctgt(i)], ... 
               [zcsrc(i),zctgt(i)]);            
    set(lin,'Color','k','LineWidth',1.2)
end
hold off, axis off vis3d


%% Plot centralities over cluster domain 

% degree, closeness and betweeness centrality
deg = M.degreeCentrality{nc};
clo = M.closenessCentrality{nc};
betm = min(M.betweenessCentrality{nc}); 
betM = max(M.betweenessCentrality{nc});

% normalized betweeness 
bet = M.betweenessCentrality{nc};
bet = (bet - betm)./(betM - betm);
bet(bet < 1e-5) = 1e-5;

% plot scale
scale1 = 800;
scale2 = 500;
scale3 = 100;

figure
subplot(3,1,1)
plotGrid(G, Ind(pucSt.(mf).compVoxelInds{nc}),...
    'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.05, ...
    'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.5);
axis off, view([82,8])

hold on
scatter3(xc,yc,zc,scale1*deg,deg,'filled');
colorbar
colormap('summer');
title('Degree centrality','fontsize',12);

%
subplot(3,1,2)
plotGrid(G, Ind(pucSt.(mf).compVoxelInds{nc}),...
    'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.05, ...
    'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.5);
axis off, view([82,8])

hold on
scatter3(xc,yc,zc,scale2*clo,clo,'filled')
colorbar
colormap('hot')
title('Closeness centrality','fontsize',12)


%
subplot(3,1,3)
plotGrid(G, Ind(pucSt.(mf).compVoxelInds{nc}),...
    'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.05, ...
    'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.5);
axis off, view([82,8])

hold on
scatter3(xc,yc,zc,scale3*bet,bet,'filled')
colormap('winter')
colorbar
title('Normalized betweeness centrality','fontsize',12)
