%% EXAMPLE: Identification of clusters in reduced domain
%
%  REMARK: This example is the head to use the MRST module co2lab
%          for CO2 injection case studies when reduced grids
%          are mandatory due to cell removal.
%           
%          In particular, we test the UNISIM-I-D model. Since it has
%          a large discontinuity interlayer region, the co2lab will 
%          remove the portion unconnected from the top, thus making
%          a reduced grid.
%           
%          That is to say, originally the grid has 93960 cells. 
%          After processGRDECL, G.cells.num returns 38466.
%          Next, after [Gt,Gaux] = topSurfaceGrid(G), 
%          Gaux.cells.num falls down to 8083.
%
%          To work with this new "reduced grid" and operate with the
%          clustering algorithm, we need to overwrite the ACTNUM
%          vector to restrain the PUC classes within the smaller domain.
%
%       
%          See "computedPUCReduced" and modifications in "maskVoxel".
%           
 

mrstVerbose off
mrstModule add mrst-gui co2lab coarsegrid matlab_bgl

case_name = 'ex6';

%% Mounting

d = DirManager(); 

%% Reservoir model

f = fullfile(d.getBenchMarksDir,'unisim-I-D','eclipse','UNISIM_I_D_ECLIPSE.DATA');
  
grdecl = readGRDECL(f);

G = processGRDECL(grdecl);
G = computeGeometry(G);

%% Rock model

% porosity
rock.poro = grdecl.PORO;
rock.poro(find(rock.poro <= 0)) = 1e-6;

% converted permeability 
rock.perm(:,1) = grdecl.PERMX * milli*darcy;
rock.perm(:,2) = grdecl.PERMY * milli*darcy;
rock.perm(:,3) = grdecl.PERMZ * milli*darcy;


%% Top surface

[Gt,Gaux] = topSurfaceGrid(G);
rock2D = averageRock(rock, Gt);
Gt.faces.neighbors = int32(Gt.faces.neighbors);
Gt.faces.nodes = int32(Gt.faces.nodes);

%% Overwrite ACTNUM

%{
Here, we will exclude the cells that do not belong to Gaux
to have a reduced grid. We will overwrite ACTNUM.

%}

grdecl2 = readGRDECL(f);
grdecl2.ACTNUM = grdecl2.ACTNUM*0;
grdecl2.ACTNUM(Gaux.cells.indexMap) = 1; % reduced grid
grdecl2.ACTNUM = logical(grdecl2.ACTNUM);
active_new = grdecl2.ACTNUM;
%save('active-co2.mat','active_new')

%% Trap analysis

% interactiveTrapping(Gt)
%trapa = trapAnalysis(Gt,fa);
%trap_volume = volumesOfTraps(Gt,trapa,unique(trapa.traps(trapa.traps>0)));

%% PUCs

% We use a fake J functional just for testing purposes
J = rock.poro.*harmmean(rock.perm,2);
J = (J - min(J))./(max(J) - min(J));

% reshape to have 3D matrix
J = reshape(J,G.cartDims);

% uses a specific version of computePUC
[PUC,nclasses,delta,divs] = computePUCReduced(J,active_new,'auto');

nofsc = 30;
puc = unique(PUC);
puc = puc(puc>0);

% demonstrative connections
C = findConnectionsByPUC(d, puc, PUC, nofsc, 'y', 1);


%% Mapping 

%{
    Here, we check if the connected components obtained over the large 
    grid domain are totally immersed into the reduced grid.
    
%}
Ind = nan(prod(G.cartDims),1);
Ind(G.cells.indexMap) = 1:G.cells.num;

P = 'PUC3';
c = 2; % choose 1,2,3, or 4 for nofsc
cvi = C.(P).compVoxelInds{c};
test = ismember(cvi,Gaux.cells.indexMap);


fprintf('IMMERSION TEST on Gaux \n')
if all(test)
    fprintf('Cluster (%s,%d) totally immersed\n',P,c);
else
    fprintf('Cluster (%s,%d) with immersion rate of: %f\n',P,c,numel(find(test))/numel(test));
    where = setdiff(cvi,Gaux.cells.indexMap);
    fprintf('Number of cells not immersed: %d\n',numel(where));a

end

%% Plot domains

%{
    Here, we make some plots just for visualization purposes.
%}

% background
plotGrid(G,'facecolor','k','facealpha',0.05,'edgecolor','none');

% top
plotGrid(Gt,'facecolor','g','facealpha',0.2,'edgecolor','k','edgealpha',0.1);

% reduced
plotGrid(Gaux,'facecolor','r','facealpha',0.2,'edgecolor','k','edgealpha',0.1);

% test cluster
% Here we could use G or Gaux, since, G contains Gaux indexing, in a
% certain sense. Ind(cvi) is necessary
plotGrid(Gaux,Ind(cvi),'facecolor','blue','edgecolor','white');
