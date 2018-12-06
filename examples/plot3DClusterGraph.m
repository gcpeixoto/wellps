%% Example: How to plot cluster centroids and its associated 3D network 
% In this tutorial, we are going to plot the cluster geometry,
% the centroid nodes and edges linking the centroids, thus
% obtaining a 3D graph representation.

%% Load grid 
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');

% compute geometry
Gc = computeGeometry(G);

%% Load connected cluster files 

% Here, we load .mat files 
load('../mat/DRT_geometric_ln_13.mat','drtSt');
load('../mat/DRT_geometric_ln_13_Metrics.mat','Maux');

% cluster ID chosen to plot 
nc = 3;

%% Inverse map to recover original indices in the correct entries
Ind = nan(prod(Gc.cartDims),1);
Ind(Gc.cells.indexMap) = 1:Gc.cells.num;

%% Plot the background grid 
% This will plot the cluster region in dimmed gray.
% You may control the colors and visual.

figure
plotGrid(Gc, Ind(drtSt.compVoxelInds{nc}),...
    'FaceColor',[0.6,0.6,0.6],'FaceAlpha',0.05, ...
    'EdgeColor',[0.6,0.6,0.6],'EdgeAlpha',0.5);
clf

%% Plot the centroids 
% To plot the cell centroids in space we need to get 
% the global indices of each cluster cell according to 'indexMap'.
% Then, we get the global stored indices from 'compVoxelInds'
% and store them as global indices in the new map.
globInd = Ind(drtSt.compVoxelInds{nc});


% centroid coordinates of cluster cells
xc = Gc.cells.centroids(globInd,1);
yc = Gc.cells.centroids(globInd,2);
zc = Gc.cells.centroids(globInd,3);
scatter3(xc,yc,zc,30,[0.6,0.2,0.1],'filled')

%% Plot the 3D graph
% First of all, we need to load the adjacency matrix of the cluster, which
% gives us the (local) indices of connected cells
adj = Maux.adjMatrix{nc};
[Iloc,Jloc] = find(adj==1);

% Once the indices are found, we create arrays to prepare the 'edge'
% plotting by choosing the source and target cells.

% source centroids
xcsrc = Gc.cells.centroids(globInd(Iloc),1); 
ycsrc = Gc.cells.centroids(globInd(Iloc),2); 
zcsrc = Gc.cells.centroids(globInd(Iloc),3); 

% target centroids
xctgt = Gc.cells.centroids(globInd(Jloc),1); 
yctgt = Gc.cells.centroids(globInd(Jloc),2); 
zctgt = Gc.cells.centroids(globInd(Jloc),3); 

hold on
% This will plot connected lines to mimic the graph edges
for i = 1:numel(Iloc)
    lin = line([xcsrc(i),xctgt(i)], ... 
               [ycsrc(i),yctgt(i)], ... 
               [zcsrc(i),zctgt(i)]);            
    set(lin,'Color','k','LineWidth',1.2)
end
hold off, axis off vis3d
