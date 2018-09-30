%% CASE: comparison of DRT clustering 

%% Load grid 
[G,PROPS] = buildModel('../benchmarks/unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE_NO_TRAILING.DATA');

% compute geometry
Gc = computeGeometry(G);

%% Parameters
% structure of several parameters (RQI, FZI, PHIZ, etc.)
P = computeParams(G,PROPS);

% field statistics for chosen properties
SA = printStats(P,{'DRTA_LN','DRTA_LOG10'},'n');
SG = printStats(P,{'DRTG_LN','DRTG_LOG10'},'n');
SN = printStats(P,{'DRTN_LN','DRTN_LOG10'},'n');

%% DRT list w/out zero

% DRT list arithmetic
drtlistA_ln = SA{1}(2:end,1);
drtlistA_log10 = SA{2}(2:end,1);

% DRT list geometric
drtlistG_ln = SG{1}(2:end,1);
drtlistG_log10 = SG{2}(2:end,1);

% DRT list normalised
drtlistN_ln = SN{1}(2:end,1);
drtlistN_log10 = SN{2}(2:end,1);

%% Mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

% cells with value
active = find(~isnan(Ind));

%% plotting DRT 3D distribution

subplot(321)
title('DRT; k-arithmetic; LN')
p1 = P.DRTA_LN(:);
p1 = p1(active);
plotCellData(G,p1,'EdgeColor','none')
axis off vis3d
colormap(jet)
colorbar

subplot(322)
title('DRT; k-arithmetic; LOG10')
p2 = P.DRTA_LOG10(:);
p2 = p2(active);
plotCellData(G,p2,'EdgeColor','none')
axis off vis3d
colormap(jet)
colorbar

subplot(323)
title('DRT; k-geometric; LN')
p3 = P.DRTG_LN(:);
p3 = p3(active);
plotCellData(G,p3,'EdgeColor','none')
axis off vis3d
colormap(jet)
colorbar

subplot(324)
title('DRT; k-geometric; LOG10')
p4 = P.DRTG_LOG10(:);
p4 = p4(active);
plotCellData(G,p4,'EdgeColor','none')
axis off vis3d
colormap(jet)
colorbar

subplot(325)
title('DRT; k-normalised; LN')
p5 = P.DRTG_LN(:);
p5 = p5(active);
plotCellData(G,p5,'EdgeColor','none')
axis off vis3d
colormap(jet)
colorbar

subplot(326)
title('DRT; k-normalised; LOG10')
p6 = P.DRTG_LOG10(:);
p6 = p6(active);
plotCellData(G,p6,'EdgeColor','none')
axis off vis3d
colormap(jet)
colorbar

%% DRT histograms 

figure 

subplot(321)
histogram(p1,'Normalization','pdf')
title('DRT; k-arithmetic; LN')

subplot(322)
histogram(p2,'Normalization','pdf')
title('DRT; k-arithmetic; LOG10')

subplot(323)
histogram(p3,'Normalization','pdf')
title('DRT; k-geometric; LN')

subplot(324)
histogram(p4,'Normalization','pdf')
title('DRT; k-geometric; LOG10')

subplot(325)
histogram(p5,'Normalization','pdf')
title('DRT; k-normalised; LN')

subplot(326)
histogram(p6,'Normalization','pdf')
title('DRT; k-normalised; LOG10')

%% DRT connections

% number of significant cells (only to save info)
nofs = 50;

drtStA_ln = findDRTConnections(drtlistA_ln, P, 'arithmetic','ln',nofs,'y', 1);
drtStA_log10 = findDRTConnections(drtlistA_log10, P, 'arithmetic','log10',nofs,'y', 1);

drtStG_ln = findDRTConnections(drtlistG_ln, P, 'geometric','ln',nofs,'y', 1);
drtStG_log10 = findDRTConnections(drtlistG_log10, P, 'geometric','log10',nofs,'y', 1);

drtStN_ln = findDRTConnections(drtlistN_ln, P, 'normalized','ln',nofs,'y', 1);
drtStN_log10 = findDRTConnections(drtlistN_log10, P, 'normalized','log10',nofs,'y', 1);

drT
%% 
