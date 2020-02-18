%% EXAMPLE 1: READ BENCHMARK GRIDS
%  Read and plot data from benchmark grids

%% Settings 

% Switch to 'on' for MRST verbose
mrstVerbose off

% Wellps classes
d = DirManager;

%% Grid reading of pre-built models

% benchmark dir
bmark = d.getBenchMarksDir;

% Choose the model
% 1: 'unisim1-eclipse'
% 2: 'unisim1-imex' (corner-point 3D view unsuported)
% 3: 'unisim2'

model = 4;
actnum = false;

switch model

    case 1
        f = fullfile(bmark,'unisim-I-D/eclipse/UNISIM_I_D_ECLIPSE.DATA');            
        actnum = true;
       
    case 2
        f = fullfile(bmark,'unisim-I-D/imex/UNISIM-I-D-POR-PERM.DAT');                       
   
    case 3
        f = fullfile(bmark,'unisim-II-D/eclipse/UNISIM_II_D_ECLIPSE.DATA');            
        actnum = true;
        
    case 4
        % nothing
        
    case 5
        % nothing
                                    
end

if model == 4
    [G,PROPS] = buildModelSPE10('mrst');
    PROPS.PHI = PROPS.PHI(:);
    PROPS.KX  = PROPS.KX(:);
else
    [G,PROPS] = buildModel(f);
end

%% Plot properties 
% 'buildModel' delivers a processed grid structure G.
% In ECLIPSE decks, MRST bypasses cells whose ACTNUM property
% is null. Hence, the number of cells in G.cells
% causes porosity and other properties in PROPS to have 
% higher dimension. To plot properties, we need to ignore
% ACTNUM. 
if actnum
    on = PROPS.ACTNUM == 1;
    PROPS.PHI = PROPS.PHI(on);
    PROPS.KX = PROPS.KX(on);
end

% plot grid  
subplot(3,1,1)
plotGrid(G)

% plot porosity
subplot(3,1,2)
plotCellData(G,PROPS.PHI);

% plot permeability
subplot(3,1,3)
plotCellData(G,PROPS.KX);


