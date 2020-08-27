%% Productivity Units 
% Script to study 3D productivity units based on proxy 
% functions of productivity potential

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


case_name = 'unisim1-prodUnits';  

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

% active cells
on = find(PROPS.ACTNUM == true);

%% Compute parameters 
P = computeParams(G,PROPS); 

%% Mapping
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

%% RQI x SO

% This section only deals with only the local grid. 
% We use the code here to plot the fields of the proxy
% function RQI x SO for visualization.

% normalized RQI
RQI_normAll = P.RQIN./max(P.RQIN);
RQI_norm = RQI_normAll(on);
ir = find(RQI_norm > 0);
RQI_norm_pos = RQI_norm(ir);

% oil sat
SO = PROPS.SO(on);
io = find(SO > 0);
SO_pos = SO(io);

% proxy: RQInorm x oil sat
RQI_SO_norm = RQI_norm(:).*SO(:);   
in = find(RQI_SO_norm > 0);
RQI_SO_norm = RQI_SO_norm(in);

% plot proxy
%Gn = extractSubgrid(G,in);
%plotCellData(Gn,RQI_SO_norm,'FaceAlpha',1.0,'EdgeColor','none')
%colormap jet

% plot RQI_norm > 0
%figure
%Gr = extractSubgrid(G,ir);
%plotCellData(Gr,RQI_norm_pos)

% plot SO > 0
%figure
%Go = extractSubgrid(G,io);
%plotCellData(Go,SO_pos)


%% Potential ranges

% We are defining a new concept based on Productivity Units (PUs)
% that are computed on the basis of a discrete function
% that defines thresholds for the proxy function. 
% This function is called here PUC (productivity unit class).
% A priori, 4 classes are expected, but only 3 really matter.
% 
% For a given reservoir cell w and proxy function F (in this case,
% F = RQI_norm(k,phi) x SO(t)), where RQI_norm(k,phi) is the 
% normalized RQI. This way, F : w -> [0,1]
%
%          { 0, if F(w) < d*
% PUC(w) = { 1, if d1 <= F(w) <  d2
%          { 2, if d2 <= F(w) <  d3
%          { 3, if d3 <= F(w) <= 1.0
%
% where d* is the average value of F for all w and d2, d3 are
% thresholds based on fixed steps:
%
% range:  [ 0 --- d* --- d2 --- d3 --- 1 ]
% class:  |   0   |   1   |   2  |   3   |   
%
%
% Therewith, we define a PU of class X as each connected cluster 
% whose all cells have PUC(w) = X.

min = mean(RQI_SO_norm); % lim inf by mean
div = linspace(min,1,4); % divisions for 3 classes

% productivity unit class
so = reshape(PROPS.SO,G.cartDims);
RQIso = RQI_normAll.*so;

PUC = RQIso;

c0 = PUC < div(1); % class 0 (ignored)
c1 = div(1) <= PUC & PUC <  div(2); % class 1
c2 = div(2) <= PUC & PUC <  div(3); % class 2
c3 = div(3) <= PUC & PUC <= div(4); % class 3

PUC(c0) = 0;
PUC(c1) = 1;
PUC(c2) = 2;
PUC(c3) = 3;

%% Productivity Units 

% compute PUs by PUC
nofsc = 10; % only for .csv
puc = 1:3; % only 3 classes matter
pucSt = findConnectionsByPUC(d,puc,PUC,nofsc,'y',1);

% PU clusters of class 3
nofc = 10; % minimum number of cells to consider to get clusters
sc1 = find(cell2mat(pucSt.PUC1.compNNodes) > nofc);
sc2 = find(cell2mat(pucSt.PUC2.compNNodes) > nofc);
sc3 = find(cell2mat(pucSt.PUC3.compNNodes) > nofc);

%% Plots

nc1 = length(sc1);
nc2 = length(sc2);
nc3 = length(sc3);
cells_sc1 = cell(1,nc1);
cells_sc2 = cell(1,nc2);
cells_sc3 = cell(1,nc3);

% plot proxy
subplot(2,2,1)
Gn = extractSubgrid(G,in);
plotCellData(Gn,RQI_SO_norm,'FaceAlpha',1.0,'EdgeColor','none')
axis off
colormap jet
title('RQI_n x so')

% PUC 1
subplot(2,2,2)
colormap jet
plotGrid(G,'FaceColor','k','FaceAlpha',0.05,'EdgeColor','none')
axis off
colorbar
title('PUC 1')
for i = 1:nc1
    cid = sc1(i);
    glob = pucSt.PUC1.compVoxelInds{cid}; % global indices
    cvim = Ind(glob); % mapped indices            
    cells_sc1{i} = cvim;
    Gg = extractSubgrid(G,cvim); % get correct indices
    plotCellData(Gg,RQIso(glob),'EdgeColor','k'); % proxy values in relaation to global
end


% PUC 2
subplot(2,2,3)
colormap jet
plotGrid(G,'FaceColor','k','FaceAlpha',0.05,'EdgeColor','none')
axis off
colorbar
title('PUC 2')
for i = 1:nc2
    cid = sc2(i);
    glob = pucSt.PUC2.compVoxelInds{cid}; % global indices
    cvim = Ind(glob); % mapped indices            
    cells_sc2{i} = cvim;
    Gg = extractSubgrid(G,cvim); % get correct indices
    plotCellData(Gg,RQIso(glob),'EdgeColor','k'); % proxy values in relaation to global
end

% PUC 3
subplot(2,2,4)
colormap jet
plotGrid(G,'FaceColor','k','FaceAlpha',0.05,'EdgeColor','none')
axis off 
colorbar
title('PUC 3')
for i = 1:nc3
    cid = sc3(i);
    glob = pucSt.PUC3.compVoxelInds{cid}; % global indices
    cvim = Ind(glob); % mapped indices            
    cells_sc3{i} = cvim;
    Gg = extractSubgrid(G,cvim); % get correct indices
    plotCellData(Gg,RQIso(glob),'EdgeColor','k'); % proxy values in relaation to global
end
