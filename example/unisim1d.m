%% TUTORIAL 

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
plotCellData(G,PROPS.PHI(on));

% plot permeability x
figure
plotCellData(G,PROPS.KX(on));

%% Compute required parameters
% Because of ACTNUM, we need to work with all the cells, thus 
% computing the required parameters even in the points where ACTNUM = 0. 
% Note that the total number of cells is the product of the 
% princiapal discretization 
ncells = prod(G.cartDims) == numel(PROPS.ACTNUM); % returns 'true'

% structure of several parameters (RQI, FZI, PHIZ, etc.)
P = computeParams(G,PROPS);
