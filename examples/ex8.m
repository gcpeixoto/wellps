%% EXAMPLE: applications with co2lab
%
% - Distance-to-trap weighting functions
% - Trap analysis and trap boundaries

mrstVerbose off
mrstModule add mrst-gui co2lab

case_name = 'ex8';

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

% COMMENTED FOR INTERACTIVE TRAP ANALYSIS TO WORK
%Gt.faces.neighbors = int32(Gt.faces.neighbors);
%Gt.faces.nodes = int32(Gt.faces.nodes);

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

%interactiveTrapping(Gt)
trapa = trapAnalysis(Gt,false);
trap_volume = volumesOfTraps(Gt,trapa,unique(trapa.traps(trapa.traps>0)));

%% Distance to trap's top

trap_i = 8; 
trap_i_top = trapa.top(10); % ideal is trapa.top(trap_i), but it seems the indices do not match in order
trap_i_cells = find(trapa.trap_regions == trap_i);
trap_i_top_centroids = Gaux.cells.centroids(trap_i_top,:);

dist_to_top = sqrt(sum(...
    [Gaux.cells.centroids(:,1) - trap_i_top_centroids(1), ...
     Gaux.cells.centroids(:,2) - trap_i_top_centroids(2), ...
     Gaux.cells.centroids(:,3) - trap_i_top_centroids(3)].^2,2));


%% Trap boundary cells

% get trap boundary cells and trap region cells
trap_i = 1;
[tcb,tcr] = getTrapBoundaryCells(Gt,trapa,trap_i);

figure
plotGrid(G,'facecolor','k','facealpha',0.1,'edgecolor','none')
plotGrid(Gt,tcr,'facecolor','g')
plotGrid(Gt,tcb,'faceColor','r')


%% Weighting functions

% --- MODEL 1: log
eta1 = 1./log(dist_to_top);


% --- MODEL 2: logistic function type 1

a = 1; b = 1; % adjust constants
D = max(dist_to_top); % D actually must be max { d(x_T,boundary) }
eta2 = exp( - a * ( dist_to_top./D ) .^ b ); 


% MODEL 3: logistic function type 2 - real logit

a = 1; b = -1e-3; c = 0; d = 0.5; % adjust constants :: change b
D = max(dist_to_top); % D actually must be max { d(x_T,boundary) }
eta3 = a./(1 + exp( - b * ( dist_to_top + c ) ) ) + d; 

%%

figure
subplot(2,3,1)
plotGrid(Gaux,'edgecolor','none','facecolor','k','facealpha',0.1)
plotGrid(Gt,trap_i_cells)
plotGrid(Gt,trap_i_top,'facecolor','red')
title('grid + trap region')


subplot(2,3,2)
plotGrid(Gt,trap_i_cells)
plotGrid(Gt,trap_i_top,'facecolor','red')
title('trap region + trap top')


subplot(2,3,3)
plotGrid(Gt,trap_i_top,'facecolor','red')
plotCellData(Gaux,eta1,'facealpha',0.5)
colorbar
title('log(d)')


subplot(2,3,4)
plotGrid(Gt,trap_i_top,'facecolor','red')
plotCellData(Gaux,eta2,'facealpha',0.5)
colorbar
title('logit_1(d)')


subplot(2,3,5)
plotGrid(Gt,trap_i_top,'facecolor','red')
plotCellData(Gaux,eta3,'facealpha',0.5)
colorbar
title('logit_pure(d)')


