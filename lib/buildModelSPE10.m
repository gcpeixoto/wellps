function [G,PROPS] = buildModelSPE10(orig)
% Read data file and build reservoir grid
%
% SYNOPSIS: 
%   [G,PROPS] = buildModelSPE10;
%
% PARAMETERS: 
%  orig  - flag used to load different data. 
%
% RETURNS:
%  G     - structure already processed in MRST for SPE10 model
%
%  PROPS - structure containing petrophysical parameters 
%        (porosity, permeability x,y,z, etc.)
%

d = DirManager;
bmark = d.getBenchMarksDir;

f = fullfile(bmark,'spe10');


% load original G structure
load(fullfile(f,'mrst/G-SPE10.mat'),'G');

% load rock structure 

switch orig % original data is in [m^2]
      
    case 'mrst'
        
        load(fullfile(f,'mrst/rock-SPE10.mat'));        
        PROPS.PHI = reshape(rock.poro,G.cartDims);
        PROPS.KX = reshape(rock.perm(:,1),G.cartDims);
        PROPS.KY = reshape(rock.perm(:,2),G.cartDims);
        PROPS.KZ = reshape(rock.perm(:,3),G.cartDims);
        
    case 'original'
        
        load(fullfile(f,'original/PHI.mat'));        
        load(fullfile(f,'original/KX.mat'));        
        load(fullfile(f,'original/KY.mat'));        
        load(fullfile(f,'original/KZ.mat'));                
        PROPS.PHI = PHI;
        PROPS.KX = KX;
        PROPS.KY = KY;
        PROPS.KZ = KZ;

end
