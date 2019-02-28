function [G,PROPS] = buildModelSPE10
% Read data file and build reservoir grid
%
% SYNOPSIS: 
%   [G,PROPS] = buildModelSPE10;
%
%
% RETURNS:
%  G     - structure already processed in MRST for SPE10 model
%
%  PROPS - structure containing petrophysical parameters 
%        (porosity, permeability x,y,z, etc.)
%


% load original G structure
load('../benchmarks/spe10/mrst/G-SPE10.mat','G');

% load rock structure 
% original data is in [m^2]
load('../benchmarks/spe10/mrst/rock-SPE10.mat');
PROPS.PHI = reshape(rock.poro,G.cartDims);
PROPS.KX = reshape(rock.perm(:,1),G.cartDims);
PROPS.KY = reshape(rock.perm(:,2),G.cartDims);
PROPS.KZ = reshape(rock.perm(:,3),G.cartDims);


