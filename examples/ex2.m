%% EXAMPLE 2: CLUSTER CENTROIDS AND 3D "SKELETON"
% In this tutorial, we are going to plot the cluster geometry,
% the centroid nodes and edges linking the centroids, thus
% obtaining a 3D graph representation.

mrstVerbose off  % turn on verbose

%% Mounting 

% class instantiation 
d = DirManager(); 

%% Grid reading
f = fullfile(d.getBenchMarksDir,'unisim-I-D','eclipse','UNISIM_I_D_ECLIPSE.DATA');

[G,PROPS] = buildModel(f);

% compute geometry
G = computeGeometry(G);

%% Load sample cluster files 

sdir = fullfile(d.getRootDir,'examples','sample');

% Here, we load the sample .mat files 
load(fullfile(sdir,'C.mat'),'drtSt');
load(fullfile(sdir,'CLinRegr.mat'),'L');
load(fullfile(sdir,'CMetrics.mat'),'M');

%% cluster ID chosen to plot 
% The maximum value is 'drtSt.allNComps'
nc = 4;

%% Inverse map to recover current logical indices at the original (unprocessed) grid
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

%% Plot the background grid 
% This will plot the cluster region in dimmed gray.
% You may control the colors and visual.

figure
plotGrid(G, Ind(drtSt.compVoxelInds{nc}),...
    'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.05, ...
    'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.5);
axis off, view([82,8])
clf

%% Plot the centroids 
% To plot the cell centroids in space we need to get 
% the global indices of each cluster cell according to 'indexMap'.
% Then, we get the global stored indices from 'compVoxelInds'
% and store them as global indices in the new map.
globInd = Ind(drtSt.compVoxelInds{nc});


% centroid coordinates of cluster cells
xc = G.cells.centroids(globInd,1);
yc = G.cells.centroids(globInd,2);
zc = G.cells.centroids(globInd,3);
scatter3(xc,yc,zc,30,[0.6,0.2,0.1],'filled')

%% Plot the 3D graph
% First of all, we need to load the adjacency matrix of the cluster, which
% gives us the (local) indices of connected cells
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
bet = M.betweenessCentrality{nc}/max(M.betweenessCentrality{nc});
bet(bet <= 0) = 1e-5;

figure
subplot(3,1,1)
plotGrid(G, Ind(drtSt.compVoxelInds{nc}),...
    'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.05, ...
    'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.5);
axis off, view([82,8])

hold on
scatter3(xc,yc,zc,10*deg,deg,'filled');
colorbar
colormap('summer');
title('Degree centrality','fontsize',12);

%
subplot(3,1,2)
plotGrid(G, Ind(drtSt.compVoxelInds{nc}),...
    'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.05, ...
    'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.5);
axis off, view([82,8])

hold on
scatter3(xc,yc,zc,100*clo,clo,'filled')
colorbar
colormap('hot')
title('Closeness centrality','fontsize',12)


%
subplot(3,1,3)
plotGrid(G, Ind(drtSt.compVoxelInds{nc}),...
    'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.05, ...
    'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.5);
axis off, view([82,8])

hold on
scatter3(xc,yc,zc,100*bet,bet,'filled')
colormap('winter')
colorbar
title('Normalized betweeness centrality','fontsize',12)



